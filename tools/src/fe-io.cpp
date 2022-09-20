
#include <config.h>

#include <stdio.h>
#include <string.h>
#include <fstream>

#include "tools-structures.h"
#include "constants.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "netcdf-proto.h"
#include "poc-time.h"
#include "rutin.h"
#include "geo.h"
#include "swap.h"

//extern double *fct_fromhms(int d, int m, double s);



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readmesh_SCHISM (const char *filename,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int i,k,n;
  int ndum,nitems,ntriangles;
  int type,ntags,tag, size;
  char keyword[1024],version[1024];
  char *line=new char[1024], *s;

  in = fopen(filename, "r");
  if(in==0) {
    return(-1);
    }

  s=fgets(line,1024,in);
  nitems=sscanf(line, "%s ", keyword);
  s=fgets(line,1024,in);
  nitems=sscanf(line, "%d %d", &mesh->ntriangles, &mesh->nvtxs);

  mesh->vertices = new vertex_t[mesh->nvtxs];

  for (i=0; i<mesh->nvtxs; i++) {
    mesh->vertices[i].code=0;
    nitems=fscanf(in, "%d %lf %lf %lf",&ndum, &mesh->vertices[i].lon, &mesh->vertices[i].lat, &mesh->vertices[i].h);
// /*------------------------------------------------------------------------------
//     convert depths into bathymetry */  
//     mesh->vertices[i].h*=-1;
    if(nitems != 4) {
      printf("error in fe_readmesh: node=%d nitems=%d\n", i,nitems);
      return(-1);
      }
    }

  mesh->triangles = new triangle_t[mesh->ntriangles];

  ntriangles=0;
  for (i=0; i<mesh->ntriangles; i++) {
    nitems=fscanf(in, "%d %d", &ndum,&type);
    if(nitems != 2) {
      printf("error in fe_readmesh\n");
      return(-1);
      }

    switch (type) {
      case 3:
/* *----------------------------------------------------------------------------
        triangle element*/
        for (k=0; k<3; k++) {
          nitems=fscanf(in, "%d", &(mesh->triangles[ntriangles].vertex[k]));
          mesh->triangles[ntriangles].vertex[k]-=1;                             //********* -1 : FORTRAN TO C index TABLE *************
          }
        ntriangles++;
        break;
      }
    }
  nitems=fscanf(in, "%s", keyword);

  mesh->ntriangles=ntriangles;
  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readmesh_GMSH_WW (const char *filename,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int i,k,n;
  int ndum,nitems,ntriangles;
  int type,ntags,tag, size;
  char keyword[1024],version[1024];

  in = fopen(filename, "r");
  if(in==0) {
    return(-1);
    }

  nitems=fscanf(in, "%s", keyword);
  nitems=fscanf(in, "%s %d %d", version, &type, &size);
  nitems=fscanf(in, "%s", keyword);

//  do { c=fgetc(in);   }  while ((c != '\n') && !feof(in));*/

  nitems=fscanf(in, "%s ", keyword);
  nitems=fscanf(in, "%d ", &mesh->nvtxs);
//  printf("#Number of nodes: %d (%d)\n",mesh->nvtxs,nitems);

  mesh->vertices = new vertex_t[mesh->nvtxs];

  for (i=0; i<mesh->nvtxs; i++) {
    mesh->vertices[i].code=0;
    nitems=fscanf(in, "%d %lf %lf %lf",&ndum, &mesh->vertices[i].lon, &mesh->vertices[i].lat, &mesh->vertices[i].h);
    if(nitems != 4) {
      printf("error in fe_readmesh: node=%d nitems=%d\n", i,nitems);
      return(-1);
      }
    }

  nitems=fscanf(in, "%s", keyword);

  nitems=fscanf(in, "%s", keyword);
  nitems=fscanf(in, "%d", &mesh->ntriangles);

  mesh->triangles = new triangle_t[mesh->ntriangles];

  ntriangles=0;
  for (i=0; i<mesh->ntriangles; i++) {
    nitems=fscanf(in, "%d %d %d", &ndum,&type,&ntags);
    if(nitems != 3) {
      printf("error in fe_readmesh\n");
      return(-1);
      }

    switch (type) {
      case 2:
/* *----------------------------------------------------------------------------
        triangle element*/

        for (k=0; k<ntags; k++) {
          nitems=fscanf(in, "%d", &tag);
          }
        for (k=0; k<3; k++) {
          nitems=fscanf(in, "%d", &(mesh->triangles[ntriangles].vertex[k]));
          mesh->triangles[ntriangles].vertex[k]-=1;                             //********* -1 : FORTRAN TO C index TABLE *************
          }
        ntriangles++;
        break;
      case 15:
/* *----------------------------------------------------------------------------
        single point element*/

        for (k=0; k<ntags; k++) {
          nitems=fscanf(in, "%d", &tag);
          }
          
        for (k=0; k<1; k++) {
          nitems=fscanf(in, "%d", &n);
          }
        mesh->vertices[n-1].code=tag+1;
        break;
      }
    }
  nitems=fscanf(in, "%s", keyword);

  mesh->ntriangles=ntriangles;
  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readmesh_GMSH (const char *filename,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  int i,k;
  int ndum,nitems,ntriangles;
  int type,ntags,tag, size;
  char dum[1024],keyword[1024],version[1024];

  in = fopen(filename, "r");
  if(in==0) {
    return(-1);
    }

  nitems=fscanf(in, "%s", keyword);
  nitems=fscanf(in, "%s %d %d", version, &type, &size);
  nitems=fscanf(in, "%s", keyword);

//  do { c=fgetc(in);   }  while ((c != '\n') && !feof(in));*/

  nitems=fscanf(in, "%s ", keyword);
  nitems=fscanf(in, "%d ", &mesh->nvtxs);
//  printf("#Number of nodes: %d (%d)\n",mesh->nvtxs,nitems);

  mesh->vertices = new vertex_t[mesh->nvtxs];

  for (i=0; i<mesh->nvtxs; i++) {
    nitems=fscanf(in, "%d %lf %lf %lf",&ndum, &mesh->vertices[i].lon, &mesh->vertices[i].lat, &mesh->vertices[i].h);
    if(nitems != 4) {
      printf("error in fe_readmesh: node=%d nitems=%d\n", i,nitems);
      return(-1);
      }
    }

  nitems=fscanf(in, "%s", keyword);

  nitems=fscanf(in, "%s", keyword);
  nitems=fscanf(in, "%d", &mesh->ntriangles);

  mesh->triangles = new triangle_t[mesh->ntriangles];

  ntriangles=0;
  for (i=0; i<mesh->ntriangles; i++) {
    nitems=fscanf(in, "%d %d %d", &ndum,&type,&ntags);
    if(nitems != 3) {
      printf("error in fe_readmesh\n");
      return(-1);
      }
    for (k=0; k<ntags; k++) {
      nitems=fscanf(in, "%d", &tag);
      }
    switch (type) {
      case 2:
/* *----------------------------------------------------------------------------
        triangle element*/
        for (k=0; k<3; k++) {
          nitems=fscanf(in, "%d", &(mesh->triangles[ntriangles].vertex[k]));
          mesh->triangles[ntriangles].vertex[k]-=1;                             //********* -1 : FORTRAN TO C index TABLE *************
          }
        ntriangles++;
        break;
      case 15:
/* *----------------------------------------------------------------------------
        single point element*/
        for (k=0; k<1; k++) {
          nitems=fscanf(in, "%d", &dum);
          }
        break;
      }
    }
  nitems=fscanf(in, "%s", keyword);

  mesh->ntriangles=ntriangles;
  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_parsemesh_GMSH (const char *filename,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------
1 2-node line.
2 3-node triangle.
3 4-node quadrangle.
4 4-node tetrahedron.
5 8-node hexahedron.
6 6-node prism.
7 5-node pyramid.
8 3-node second order line (2 nodes associated with the vertices and 1 with the edge).
9 6-node second order triangle (3 nodes associated with the vertices and 3 with the edges).
10 9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face).
11 10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges).
12 27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume).
13 18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces).
14 14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face).
15 1-node point.
16 8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges).
17 20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges).
18 15-node second order prism (6 nodes associated with the vertices and 9 with the edges).
19 13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges).
20 9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)
21 10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)
22 12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)
23 15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)
24 15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)
25 21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)
26 4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)
27 5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)
28 6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)
29 20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)
30 35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)
31 56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)
-----------------------------------------------------------------------------*/
{
  FILE *in=NULL;
  double x,y,z;
  int i,k;
  int ndum,nitems;
  int ntriangles;
  int type,ntags,tag, size;
  char dum[1024],keyword[1024],version[1024];

  string gmsh_type[50];
  gmsh_type[1] =string("2-node line");
  gmsh_type[2] =string("3-node triangle");
  gmsh_type[3] =string("4-node quadrangle");
  gmsh_type[4] =string("4-node tetrahedron");
  gmsh_type[5] =string("8-node hexahedron");
  gmsh_type[6] =string("6-node prism");
  gmsh_type[7] =string("5-node pyramid");
  gmsh_type[8] =string("3-node second order line (2 nodes associated with the vertices and 1 with the edge)");
  gmsh_type[9] =string("6-node second order triangle (3 nodes associated with the vertices and 3 with the edges)");
  gmsh_type[10]=string("9-node second order quadrangle (4 nodes associated with the vertices, 4 with the edges and 1 with the face)");
  gmsh_type[11]=string("10-node second order tetrahedron (4 nodes associated with the vertices and 6 with the edges)");
  gmsh_type[12]=string("27-node second order hexahedron (8 nodes associated with the vertices, 12 with the edges, 6 with the faces and 1 with the volume)");
  gmsh_type[13]=string("18-node second order prism (6 nodes associated with the vertices, 9 with the edges and 3 with the quadrangular faces)");
  gmsh_type[14]=string("14-node second order pyramid (5 nodes associated with the vertices, 8 with the edges and 1 with the quadrangular face)");
  gmsh_type[15]=string("1-node point");
  gmsh_type[16]=string("8-node second order quadrangle (4 nodes associated with the vertices and 4 with the edges)");
  gmsh_type[17]=string("20-node second order hexahedron (8 nodes associated with the vertices and 12 with the edges)");
  gmsh_type[18]=string("15-node second order prism (6 nodes associated with the vertices and 9 with the edges)");
  gmsh_type[19]=string("13-node second order pyramid (5 nodes associated with the vertices and 8 with the edges)");
  gmsh_type[20]=string("9-node third order incomplete triangle (3 nodes associated with the vertices, 6 with the edges)");
  gmsh_type[21]=string("10-node third order triangle (3 nodes associated with the vertices, 6 with the edges, 1 with the face)");
  gmsh_type[22]=string("12-node fourth order incomplete triangle (3 nodes associated with the vertices, 9 with the edges)");
  gmsh_type[23]=string("15-node fourth order triangle (3 nodes associated with the vertices, 9 with the edges, 3 with the face)");
  gmsh_type[24]=string("15-node fifth order incomplete triangle (3 nodes associated with the vertices, 12 with the edges)");
  gmsh_type[25]=string("21-node fifth order complete triangle (3 nodes associated with the vertices, 12 with the edges, 6 with the face)");
  gmsh_type[26]=string("4-node third order edge (2 nodes associated with the vertices, 2 internal to the edge)");
  gmsh_type[27]=string("5-node fourth order edge (2 nodes associated with the vertices, 3 internal to the edge)");
  gmsh_type[28]=string("6-node fifth order edge (2 nodes associated with the vertices, 4 internal to the edge)");
  gmsh_type[29]=string("20-node third order tetrahedron (4 nodes associated with the vertices, 12 with the edges, 4 with the faces)");
  gmsh_type[30]=string("35-node fourth order tetrahedron (4 nodes associated with the vertices, 18 with the edges, 12 with the faces, 1 in the volume)");
  gmsh_type[31]=string("56-node fifth order tetrahedron (4 nodes associated with the vertices, 24 with the edges, 24 with the faces, 4 in the volume)");

  in = fopen(filename, "r");
  if(in==0) {
    return(-1);
    }

/* *----------------------------------------------------------------------------
  mesh format*/
  nitems=fscanf(in, "%s", keyword);
  nitems=fscanf(in, "%s %d %d", version, &type, &size);
  nitems=fscanf(in, "%s", keyword);

/* *----------------------------------------------------------------------------
  nodes*/
  nitems=fscanf(in, "%s ", keyword);
  nitems=fscanf(in, "%d ", &mesh->nvtxs);

  for (i=0; i<mesh->nvtxs; i++) {
    nitems=fscanf(in, "%d %lf %lf %lf",&ndum, &x,&y,&z);
    if(nitems != 4) {
      printf("error in fe_readmesh: node=%d nitems=%d\n", i,nitems);
      return(-1);
      }
    }
  nitems=fscanf(in, "%s", keyword);

/* *----------------------------------------------------------------------------
  elements*/
  nitems=fscanf(in, "%s", keyword);
  nitems=fscanf(in, "%d", &mesh->ntriangles);

  ntriangles=0;
  for (i=0; i<mesh->ntriangles; i++) {
    nitems=fscanf(in, "%d %d %d", &ndum,&type,&ntags);
    if(nitems != 3) {
      printf("error in fe_readmesh\n");
      return(-1);
      }
    for (k=0; k<ntags; k++) {
      nitems=fscanf(in, "%d", &tag);
      }
    switch (type) {
      case 2:
/* *----------------------------------------------------------------------------
        triangle element*/
        for (k=0; k<3; k++) {
          nitems=fscanf(in, "%d", &dum);
          }
        ntriangles++;
        break;
      case 15:
/* *----------------------------------------------------------------------------
        single point element*/
        for (k=0; k<1; k++) {
          nitems=fscanf(in, "%d", &dum);
          }
        break;
      }
    }
  nitems=fscanf(in, "%s", keyword);

  mesh->ntriangles=ntriangles;
  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_readmesh_QUODDY(const char *ElementFile, const char *NodeFile, const char *DepthFile,mesh_t *mesh, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in1=NULL,*in2=NULL;
  char *line=new char[256],*s __attribute__((unused));
  int i,m,n;
  int nitems,n1,n2,n3,count,status,smallest;
  char numbering;

/*------------------------------------------------------------------------------
  element and node list format*/
  in1 = fopen(ElementFile, "r");
  if(in1==0)
      return(-1);

  in2 = fopen(NodeFile, "r");
  if(in2==0)
      return(-1);

/*------------------------------------------------------------------------------
  rough estimate of the maximum count of nodes*/
  n=0;
  count=0;
  smallest=1;
  while(!feof(in1)) {
    nitems=fscanf(in1, "%d %d %d", &n1,&n2,&n3);
    if(nitems !=3) break;
    count++;
    updatemax(&n,n1);
    updatemax(&n,n2);
    updatemax(&n,n3);
    updatemin(&smallest,n1);
    updatemin(&smallest,n2);
    updatemin(&smallest,n3);
    }

  if(smallest==0) {
    mesh->nvtxs=n+1;
    numbering='C';
    }
  else {
    mesh->nvtxs=n;
    numbering='F';
    }
    
  mesh->ntriangles=count;
  if(debug)  {
    printf("#Number of elements: %6d \n",mesh->ntriangles);
    printf("#Number of nodes   : %6d (numbering=%c)\n",mesh->nvtxs,numbering);
    }

  exitIfNull(mesh->triangles = new triangle_t[mesh->ntriangles]);

  rewind(in1);
  for(m=0;m<mesh->ntriangles;m++) {
    nitems=fscanf(in1, "%d %d %d", &n1,&n2,&n3);
    if(nitems !=3) break;
    if(smallest==1) { n1-- ; n2-- ; n3-- ; }
    mesh->triangles[m].vertex[0]=n1;
    mesh->triangles[m].vertex[1]=n2;
    mesh->triangles[m].vertex[2]=n3;
    }

  exitIfNull(mesh->vertices=new vertex_t[mesh->nvtxs]);

  status=fe_e2n(mesh);

  for (i=0; i<mesh->nvtxs; i++) {
    s=fgets(line,256,in2);
    int dum;
/*------------------------------------------------------------------------------
    warning: '\n' will be accounted as token */
    nitems=count_token(line);
    if(nitems==3) {
      nitems=sscanf(line, "%lf %lf", &mesh->vertices[i].lon,&mesh->vertices[i].lat);
      }
    else if(nitems==4) {
      nitems=sscanf(line, "%d %lf %lf", &dum, &mesh->vertices[i].lon,&mesh->vertices[i].lat);
      }
    else TRAP_ERR_EXIT(status,"%s: error reading %s",__func__,NodeFile);
//     do {
//       c=fgetc(in2);
//       }  while ((c != '\n') && !feof(in2));
    }

  fclose(in1);
  fclose(in2);
  
  if(DepthFile!=0) {
    in2 = fopen(DepthFile, "r");
    if(in2==0)
      return(-1);
    for (i=0; i<mesh->nvtxs; i++) {
      s=fgets(line,256,in2);
      double x,y;
      float z;
      nitems=sscanf(line, "%lf %lf %f",&x, &y, &z);
      if(nitems<3) nitems=sscanf(line, "%f",&z);
      mesh->vertices[i].h=z;
      }
    fclose(in2);
    }
  
  status=fe_initaffine(mesh);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readmesh_TGL (const char *ElementFile, const char *NodeFile, mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  FILE *in1=NULL,*in2=NULL;
  int i,m,n,nelts,nvtxs,dum;
  int nitems,n1,n2,n3,count,row,status;
  char c;

/*------------------------------------------------------------------------------
  element and node list format*/

  in1 = fopen(ElementFile, "r");
  if(in1==0)
      return(-1);

  in2 = fopen(NodeFile, "r");
  if(in2==0)
      return(-1);

/*------------------------------------------------------------------------------
  rough estimate of the maximum count of nodes*/
  nitems=fscanf(in1, "%d %d %d", &nelts,&nvtxs,&dum);
  n=0;
  count=0;
  while(!feof(in1)) {
    nitems=fscanf(in1, "%d %d %d %d", &row, &n1,&n2,&n3);
    if(nitems !=4) break;
    count++;
    updatemax(&n,n1);
    updatemax(&n,n2);
    updatemax(&n,n3);
    }

  mesh->nvtxs=n+1;
  mesh->ntriangles=count;
//  printf("#Number of nodes: %d \n",mesh->nvtxs);

  mesh->triangles=(triangle_t *) new triangle_t[mesh->ntriangles];

  rewind(in1);
  nitems=fscanf(in1, "%d %d %d", &nelts,&nvtxs,&dum);
  for(m=0;m<mesh->ntriangles;m++) {
    nitems=fscanf(in1, "%d %d %d %d", &row, &n1,&n2,&n3);
    if(nitems !=4) break;
    mesh->triangles[m].vertex[0]=n1;
    mesh->triangles[m].vertex[1]=n2;
    mesh->triangles[m].vertex[2]=n3;
    }

  mesh->vertices= new vertex_t[mesh->nvtxs];
  if(mesh->vertices == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }

/* *-----------------------------------------------------------------------------
  create neighbour array from element array */
  for (i=0; i<mesh->nvtxs; i++) {
    mesh->vertices[i].null_value();
    }
  status=fe_e2n(mesh);

  nitems=fscanf(in2, "%d %d %d %d",&nvtxs,&dum,&dum,&dum);
  if(nitems !=4) goto error;
  for (i=0; i<mesh->nvtxs; i++) {
    nitems=fscanf(in2, "%d %lf %lf", &row, &mesh->vertices[i].lon,&mesh->vertices[i].lat);
    do { c=fgetc(in2);   }  while ((c != '\n') && !feof(in2));
    }

  fclose(in1);
  fclose(in2);
  
  status=fe_initaffine(mesh);

  return(0);

  error:
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readmesh_NGH(const char *filename,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  vertex_t *set=NULL;
  int i,j,nndes,maxn,ndum=-1;
  int *LocalNgh=NULL;
  int status=-1, srcLine, nitems, LocalNngh;
  frame_t frame;

/*------------------------------------------------------------------------------
  neighbour list format*/

  in = fopen(filename, "r");
  if(in==0) {
    TRAP_ERR_RETURN(-1,1,"cannot open %s (%d %s)\n", filename, errno, strerror(errno));
    }

  nitems=fscanf(in, "%d ", &mesh->nvtxs);
  if(nitems!=1) {
    srcLine=__LINE__;goto error;
    }
  
//  printf("#Number of nodes: %d \n",mesh->nvtxs);
  nitems=fscanf(in, "%d ", &mesh->nnghm);
  if(nitems!=1) {
    srcLine=__LINE__;goto error;
    }
    
  nitems=fscanf(in, "%lf %lf %lf %lf", &frame.xmax,&frame.ymax,&frame.xmin,&frame.ymin);
  if(nitems!=4) {
    srcLine=__LINE__;goto error;
    }

  nndes=mesh->nvtxs;
  maxn=mesh->nnghm;

  mesh->vertices= new vertex_t[nndes];
  for (i=0; i<nndes; i++) mesh->vertices[i].null_value();

  set= mesh->vertices;

  if(set == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }

  LocalNgh = new int[maxn];
  for (i=0; i<nndes; i++) {
    nitems=fscanf(in, "%d %lf %lf %d %lf", &ndum,&set[i].lon,&set[i].lat,&set[i].code,&set[i].h);
    if(nitems!=5) {
      srcLine=__LINE__;goto error;
      }
/* *------------------------------------------------------------------------
    conversion to cm, inadequate*/
//    set[i].h*=100.;
    set[i].c=cos(set[i].lat*d2r);
    LocalNngh=0;
    for(j=0; j<maxn; j++) {
      nitems=fscanf(in, "%d ", &LocalNgh[j]);
      if(nitems!=1) {
        printf("%d %lf %lf %d %lf\n", ndum,set[i].lon,set[i].lat,set[i].code,set[i].h);
        printf("reading %d th neighbour of node %d (over %d)\n",j,i,maxn);
        srcLine=__LINE__;goto error;
        }
/*------------------------------------------------------------------------------
      neighbour index correction*/
      LocalNgh[j]--;
      if(LocalNgh[j]>-1){LocalNngh=j+1;}
      }
    set[i].nngh=LocalNngh;
    if(LocalNngh<=0)
      continue;
    set[i].ngh=new int[set[i].nngh];
    for(j=0; j<set[i].nngh; j++) set[i].ngh[j]=LocalNgh[j];
    }
  
  status=0;
  
error:
  fclose(in);
  deletep(&LocalNgh);
  
  if(status){
    fprintf(stderr,"%s:%d: %d items found for %d in %s\n",strrchr0(__FILE__,'/'),srcLine,nitems,ndum,filename);
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readmesh_TELEMAC (const char *rootname, mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in1=NULL,*in2=NULL;
  int i,m,n,dum;
  int fmt_ele,fmt_nod;
  int nitems,n1,n2,n3,count,status;
  char *elements=NULL,*nodes=NULL;
  char  line[1000];

  elements=new char[strlen(rootname)+4];
  nodes=new char[strlen(rootname)+4];

  sprintf(elements,"%s.ele",rootname);
  sprintf(nodes,"%s.nod",rootname);

/*------------------------------------------------------------------------------
  element and node list format*/
  in1 = fopen(elements, "r");
  if(in1==0)
      return(-1);
  
  fgets(line,999,in1);
  nitems=sscanf(line, "%d %d %d %d", &dum, &n1, &n2, &n3);
  switch(nitems) {
    case 3:
      fmt_ele=0;
      break;
    case 4:
      fmt_ele=1;
      break;
    default:
      return(-1);
      break;
    }
  rewind (in1);
  
  in2 = fopen(nodes, "r");
  if(in2==0)
      return(-1);

  fgets(line,999,in2);
  if(strcmp(line,"NXY")==0) {
    fmt_nod=0;
    }
  else if(strncmp(line,"XYZ",3)==0) {
    fmt_nod=1;
    }
  else {
    fmt_nod=0;
    rewind (in2);
    }

/*------------------------------------------------------------------------------
  rough estimate of the maximum count of nodes*/
  n=0;
  count=0;
  while(!feof(in1)) {
    switch(fmt_ele) {
      case 0:
        nitems=fscanf(in1, "%d %d %d", &n1, &n2, &n3);
        if(nitems !=3) break;
        count++;
        break;
      case 1:
        nitems=fscanf(in1, "%d %d %d %d", &dum, &n1, &n2, &n3);
        if(nitems !=4) break;
        count++;
        break;
      }
    updatemax(&n,n1);
    updatemax(&n,n2);
    updatemax(&n,n3);
    }

  mesh->nvtxs=n;
  mesh->ntriangles=count;
//  printf("#Number of nodes: %d \n",mesh->nvtxs);

  exitIfNull(mesh->triangles = new triangle_t[mesh->ntriangles]);

  rewind(in1);
  for(m=0;m<mesh->ntriangles;m++) {
    switch(fmt_ele) {
      case 0:
        nitems=fscanf(in1, "%d %d %d", &n1, &n2, &n3);
        if(nitems !=3) break;
        count++;
        break;
      case 1:
        nitems=fscanf(in1, "%d %d %d %d", &dum, &n1, &n2, &n3);
        if(nitems !=4) break;
        count++;
        break;
      }
    mesh->triangles[m].vertex[0]=n1-1;
    mesh->triangles[m].vertex[1]=n2-1;
    mesh->triangles[m].vertex[2]=n3-1;
    }

  exitIfNull(mesh->vertices=new vertex_t[mesh->nvtxs]);

  status=fe_e2n(mesh);

  for (i=0; i<mesh->nvtxs; i++) {
    switch(fmt_nod) {
      case 0:
        nitems=fscanf(in2, "%d %lf %lf", &count, &mesh->vertices[i].lon,&mesh->vertices[i].lat);
        if(nitems !=3) break;
        break;
      case 1:
        nitems=fscanf(in2, "%lf %lf %f", &mesh->vertices[i].lon,&mesh->vertices[i].lat,&mesh->vertices[i].h);
        if(nitems !=3) break;
        break;
      }
//     do {
//       c=fgetc(in2);
//       }  while ((c != '\n') && !feof(in2));
    }

  fclose(in1);
  fclose(in2);

  return(0);
}

void CONV_LAMBERT_DEGDEC_TAB(const int NPOIN, double* X , double* Y , double* LAMBDATAB , double* PHITAB )
{

// C***********************************************************************
// C  TELEMAC 2D VERSION 6.0   25/06/10  C-T PHAM (LNHE) 01 30 87 85 93
// C
// C***********************************************************************
// C
// C FONCTION  : CONVERSION DE COORDOONNEES METRIQUES LAMBERT 1 NORD
// C             EN LATITUDES, LONGITUDES (DEGRES DECIMAUX)
// C
// C-----------------------------------------------------------------------
// C
// C FUNCTION: CONVERSION OF COORDINATES METRIC LAMBERT 1 NORTH
// C           INTO LATITUDES, LONGITUDES (DECIMAL DEGREES)
// C
// C-----------------------------------------------------------------------
// C
// C-----------------------------------------------------------------------
// C                             ARGUMENTS
// C .________________.____.______________________________________________.
// C |      NOM       |MODE|                   ROLE                       |
// C |________________|____|______________________________________________|
// C |   X,Y          |<-->| METRIC COORDONINATES (LAMBERT 1 NORTH)
// C |   LAMBDA       |<-->| LONGITUDE
// C |   PHI          |<-->| LATITUDE
// C |________________|____|______________________________________________|
// C MODE : -->(DONNEE NON MODIFIEE), <--(RESULTAT), <-->(DONNEE MODIFIEE)
// C
// C-----------------------------------------------------------------------
/* *--------------------------------------------------------------------------------------------------------------------------------------------------
 Projection  		Système  	Parallèle 1  	Parallèle 2  	Latitude origine  	Longitude origine  	X0  		Y0
 Lambert I 		NTF 	48°35'54.682" N 	50°23'45.282" N 	55 gr N 	2°20'14.025" E 		600 000 	1 200 000
 Lambert II 		NTF 	45°53'56.108" N 	47°41'45.652" N 	52 gr N 	2°20'14.025" E 		600 000 	2 200 000
 Lambert II étendu 	NTF 	45°53'56.108" N 	47°41'45.652" N 	52 gr N 	2°20'14.025" E 		600 000 	2 200 000
 Lambert III 		NTF 	43°11'57.449" N 	44°59'45.938" N 	49 gr N 	2°20'14.025" E 		600 000 	3 200 000
 Lambert IV 		NTF 	41°33'37.396" N 	42°46'03.588" N 	43.85 gr N 	2°20'14.025" E 		234.358 	4 185 861.369
 Lambert 93 		RGF93 	44°00'00.000" N 	49°00'00.000" N 	46°30' N 	3°00'00.000" E 		700 000 	6 600 000
----------------------------------------------------------------------------------------------------------------------------------------------------*/
      double pi=M_PI;
      double Z,RRR;
      double EEE,LAMBDAC,NNN,CCC,XS,YS;
      double GAMMA;
      double LAMBDA,PHI;
      double PHIM,LATISO,ES2,EPSILON;

      double AAA,GNORM;
      double HE;
      double TX,TY,TZ;

      double CORRPHI;

      int I,J;

//      PI = ACOS(-1.); /* deja ? */

      EPSILON = 1.E-11;

// LAMBERT 1 NORD

// LAMBDAC : MERIDIEN DE PARIS / MERIDIEN DE GREENWICH

      LAMBDAC = (2.+20./60.+14.025/3600.)*d2r;

// TRANSFO PROJECTION CC LAMBERT X, Y --> COORD GEO LAMBERT LAMBDA, PHI
// (IGN : ALG0004)

      EEE = 0.08248325676;
//
      ES2 = EEE/2.;
//
      NNN = 0.760405966;   // NNN = 0.7604059656
      CCC = 11603796.9767; // CCC = 11603796.98
      XS  = 600000.;
      YS  = 5657616.674;

      for (J=0;J<NPOIN;J++) {
        RRR = sqrt((X[J]-XS)*(X[J]-XS)+(Y[J]-YS)*(Y[J]-YS));
        GAMMA = atan((X[J]-XS)/(YS-Y[J]));

        LAMBDA = GAMMA/NNN+LAMBDAC;

        LATISO = -1./NNN*log(abs(RRR/CCC));

// CALCUL DE LA LAT PHI A PART DE LA LAT ISO LATISO
// (IGN ALG0002)

// I = 0
        PHIM = 2.*atan(exp(LATISO))-pi/2.;
// I = 1

       double dummy=exp(ES2*log( (1.+EEE*sin(PHIM))/(1.-EEE*sin(PHIM))));
        PHI  = 2.*atan(exp(LATISO)*dummy)-pi/2.;

      I = 1;

       while (abs(PHI-PHIM)>=EPSILON) {
          PHIM = PHI;
          dummy=exp(ES2*log( (1.+EEE*sin(PHIM))/(1.-EEE*sin(PHIM))));
          PHI  = 2.*atan(exp(LATISO)*dummy)-pi/2.;
          I = I + 1;
        }

// TRANSFO COORD GEO LAMBDA, PHI --> COORD CARTESIENNES X, Y, Z
// (IGN : ALG0009)

// CALCUL DE LA GRANDE NORMALE GNORM PREALABLE
// LAT PHI --> GNORM (IGN : ALG0021)

      HE = 0.; // POUR LE TEST, ET PAS UTILE POUR LA FINALITE

      AAA = 6378249.2;
      EEE = 0.08248325679;

      GNORM = AAA/sqrt(1.-EEE*EEE*sin(PHI)*sin(PHI));

      X[J] = (GNORM+HE)*cos(PHI)*cos(LAMBDA);
      Y[J] = (GNORM+HE)*cos(PHI)*sin(LAMBDA);
      Z    = (GNORM*(1.-EEE*EEE)+HE)*sin(PHI);

// TRANSFO COORD A 7 PARAM ENTRE 2 SYST GEODESIQUES
// TRANSFO SIMPLIFIEE, JUSTE UNE TRANSLATION A 3 PARAM (IGN :  ALG0013 SIMPLIFIE)
// NTF LAMBERT --> WGS84 (IGN)

      TX = -168.;
      TY = -60.;
      TZ =  320.;

      X[J] = X[J] + TX;
      Y[J] = Y[J] + TY;
      Z    = Z + TZ;

// TRANSFO COORD CART X, Y, Z --> COORD GEO LAMBDA, PHI
// METHODE DE HEISKANEN-MORITZ-BOUCHER (IGN : ALG0012)
 
      LAMBDA = atan(Y[J]/X[J]);
// I = 0
      PHIM   = atan(Z/(sqrt(X[J]*X[J]+Y[J]*Y[J])*(1.-(AAA*EEE*EEE)/(sqrt(X[J]*X[J]+Y[J]*Y[J]+Z*Z)))));
// I = 1
      PHI   = atan(Z/(sqrt(X[J]*X[J]+Y[J]*Y[J])*(1.-(AAA*EEE*EEE*cos(PHIM))/(sqrt(X[J]*X[J]+Y[J]*Y[J])*sqrt(1.-EEE*EEE*sin(PHIM)*sin(PHIM))))));

      I = 1;

      while (abs(PHI-PHIM)>=EPSILON) {
        PHIM = PHI;
        PHI   = atan(Z/(sqrt(X[J]*X[J]+Y[J]*Y[J])*(1.-(AAA*EEE*EEE*cos(PHIM))/(sqrt(X[J]*X[J]+Y[J]*Y[J])*sqrt(1.-EEE*EEE*sin(PHIM)*sin(PHIM))))));
        I = I + 1;
      }

// CORRECTION D ANGLE SUR PHI POUR RECALAGE

      CORRPHI = -5.439609E-5;
      PHI = PHI + CORRPHI;

// CONVERSION EN DEGRES DECIMAUX

      LAMBDA = LAMBDA*r2d;
      PHI    = PHI*r2d;

      LAMBDATAB[J] = LAMBDA;
      PHITAB[J]    = PHI;
  }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_readmesh_TELEMAC_BINARY (const char *filename, mesh_t *mesh, const char *proj4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  int i,j,k,status;
  char *title=NULL,*name=NULL,**varname=NULL;
  int recsize,nval;
  int nframes;
  int *ibuf=NULL;
  float *rbuf=NULL;
  int *nbv=NULL;
  int cartesian=1;
  bool need_swap=0;
  double t,*time;
  date_t start;
  
  in = fopen(filename, "r");
  if(in==0) {
    return(-1);
    }

/* *----------------------------------------------------------------------------
  Title */
  fread(&recsize,sizeof(int),1,in);
  
  if(recsize<0)   need_swap=1;
  if(recsize>256) need_swap=1;
  
  if(need_swap) {
    printf("TELEMAC_BINARY, binary will be swapped\n");
    }

  if(need_swap) {
    recsize=lswap(recsize);
    }

  title=new char[recsize+1];
  fread(title,recsize,1,in);
  title[recsize]=0;
  
  printf("TELEMAC_BINARY, title: %s (%d)\n",title,recsize);

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

/* *----------------------------------------------------------------------------
  NBV */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(int);
  nbv=new int[nval];
  printf("TELEMAC_BINARY, variables nrecords : %d\n",nval);

  fread(nbv,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) nbv[i]=lswap(nbv[i]);
    }
  for(i=0;i<nval;i++) printf("TELEMAC_BINARY, variable %d  : %d\n",i,nbv[i]);


  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

/* *----------------------------------------------------------------------------
  NAMES */
  varname=new char*[nbv[0]];
  for(i=0;i<nbv[0];i++) {
    fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }

    name=new char[recsize+1];
    fread(name,recsize,1,in);
    name[recsize]=0;
    
    printf("TELEMAC_BINARY, variable name %d: %s (%d)\n",i,name,recsize);
    varname[i]=strdup(name);
    
    delete[] name;

    fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }
    }

/* *----------------------------------------------------------------------------
  IPARAM */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(int);
  ibuf=new int[nval];
  printf("TELEMAC_BINARY, parameters nrecords : %d\n",nval);

  fread(ibuf,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
    }

  for(i=0;i<nval;i++) printf("TELEMAC_BINARY, parameter %d  : %d\n", i, ibuf[i]);

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  int date=(ibuf[9]==1);

  delete[] ibuf;
  
  if(date==1) {
    fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }

    nval=recsize/sizeof(int);
    ibuf=new int[nval];
    printf("TELEMAC_BINARY, date nrecords: %d\n",nval);

    fread(ibuf,recsize,1,in);
    if(need_swap) {
      for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
      }

    for(i=0;i<nval;i++) printf("TELEMAC_BINARY, date record %d  : %d\n",i,ibuf[i]);

    fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }

    start.year  =ibuf[0];
    start.month =ibuf[1];
    start.day   =ibuf[2];
    start.second=ibuf[3]*3600+ibuf[4]*60+ibuf[5];
  
    delete[] ibuf;
    }
  
/* *----------------------------------------------------------------------------
  CARDINAL */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(int);
  ibuf=new int[nval];
  printf("TELEMAC_BINARY, cardinal nrecords : %d\n",nval);

  fread(ibuf,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
    }
  for(i=0;i<nval;i++) printf("TELEMAC_BINARY, cardinal %d  : %d\n",i,ibuf[i]);

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  mesh->ntriangles=ibuf[0];
  mesh->nvtxs=ibuf[1];

  mesh->vertices = new vertex_t[mesh->nvtxs];
  mesh->triangles = new triangle_t[mesh->ntriangles];

  delete[] ibuf;

/* *----------------------------------------------------------------------------
  CONNECTIVITY */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(int);
  ibuf=new int[nval];

  fread(ibuf,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
    }

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }
  for (i=0; i<mesh->ntriangles; i++) {
    for (k=0; k<3; k++) {
      mesh->triangles[i].vertex[k]=ibuf[i*3+k]-1;
      }
    }

  delete[] ibuf;

/* *----------------------------------------------------------------------------
  BOUNDARY CODES */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(int);
  ibuf=new int[nval];

  fread(ibuf,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
    }

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }
  for(i=0;i<mesh->nvtxs;i++) {
    mesh->vertices[i].code=ibuf[i];
    }

  delete[] ibuf;

/* *----------------------------------------------------------------------------
  LONGITUDES */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(float);
  rbuf=new float[nval];

  fread(rbuf,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) rbuf[i]=fswap(rbuf[i]);
    }

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }
  for(i=0;i<mesh->nvtxs;i++) {
/* *----------------------------------------------------------------------------
    !!! translation en x */
//    mesh->vertices[i].lon=rbuf[i]-600000.;
    mesh->vertices[i].lon=rbuf[i];
    }
  delete[] rbuf;

/* *----------------------------------------------------------------------------
  LATITUDES */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(float);
  rbuf=new float[nval];

  fread(rbuf,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) rbuf[i]=fswap(rbuf[i]);
    }

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }
  for(i=0;i<mesh->nvtxs;i++) {
    mesh->vertices[i].lat=rbuf[i];
    }
  delete[] rbuf;

  float *z=new float[mesh->nvtxs];
  float *u=new float[mesh->nvtxs];
  float *v=new float[mesh->nvtxs];
  
  FILE *sample=fopen("sample.dat","w");
  
  nframes=0;
  time=new double[1000];
  double tprev=0;
  
  while(true) {
/* *----------------------------------------------------------------------------
    TIME */
    int nitems;
    nitems=fread(&recsize,sizeof(int),1,in);
    if(nitems!=1) break;
  
    if(need_swap) {
      recsize=lswap(recsize);
      }

    nval=recsize/sizeof(float);
    printf("TELEMAC_BINARY, time nrecords: %d %d\n",nval,nframes);
    if(nval==0) {
      printf("empty time record at frame %d\n",nframes);
//      break;
      recsize=4;
      nval=1;
      }
    rbuf=new float[nval];

    nitems=fread(rbuf,recsize,1,in);
    if(need_swap) {
      for(i=0;i<nval;i++) rbuf[i]=fswap(rbuf[i]);
      }
    
    double t0=cnes_time(start, 'd');
  
    for(i=0;i<nval;i++) {
      printf("TELEMAC_BINARY, time %d  : %lf\n",nframes,rbuf[i]);
      time[0]=t0+rbuf[i]/(24.*3600.);
      if((rbuf[i]==0.)&&(nframes>0)) time[0]=tprev+1800/(24.*3600.);
      tprev=time[0];
      }
    
    nitems=fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }
    delete[] rbuf;
  
/* *----------------------------------------------------------------------------
    VARIABLES */
    for(j=0;j<nbv[0];j++) {
      nitems=fread(&recsize,sizeof(int),1,in);
      if(nitems!=1) break;
      if(need_swap) {
        recsize=lswap(recsize);
        }

      nval=recsize/sizeof(float);
      printf("TELEMAC_BINARY, variable %d nrecords: %d\n",j,nval);
      if(nval==0) {
        long pos;
        pos=ftell(in);
        printf("empty record at variable %d, frame %d (pos=%ld bytes)\n",j,nframes, pos);
//        break;
        recsize=4*mesh->nvtxs;
        nval=mesh->nvtxs;
        }
      rbuf=new float[nval];

      nitems=fread(rbuf,recsize,1,in);
      if(nitems!=1) break;
      if(need_swap) {
        for(i=0;i<nval;i++) rbuf[i]=fswap(rbuf[i]);
        }

      nitems=fread(&recsize,sizeof(int),1,in);
      if(nitems!=1) break;
      if(need_swap) {
        recsize=lswap(recsize);
        }
      if(strstr(varname[j],"FOND")!=0 or strstr(varname[j],"BOTTOM")!=0) {
        for(i=0;i<mesh->nvtxs;i++) {
          mesh->vertices[i].h=rbuf[i];
          }
        }
      if(strstr(varname[j],"SURFACE LIBRE")!=0) {
        for(i=0;i<mesh->nvtxs;i++) {
          z[i]=rbuf[i];
          }
        }
      if(strstr(varname[j],"VITESSE U")!=0) {
        for(i=0;i<mesh->nvtxs;i++) {
          u[i]=rbuf[i];
          }
        }
      if(strstr(varname[j],"VITESSE V")!=0) {
        for(i=0;i<mesh->nvtxs;i++) {
          v[i]=rbuf[i];
          }
        }
      delete[] rbuf;
      }
    i=0;
    fprintf(sample,"%12.6lf %9.4f %9.4f %9.4f\n",time[0],z[i],u[i],v[i]);
    nframes++;
    }

  fclose(sample);
  
  status=fe_e2n(mesh);
  fclose(in);
  
  cartesian=0;
  for(i=0;i<mesh->nvtxs;i++) {
    t=mesh->vertices[i].lon;
    if(t>360.0) {
      cartesian=1;
      break;
      }
    if(t<-360.0) {
      cartesian=1;
      break;
      }
    }
    
  if(cartesian==1) {
//     projPJ ref=NULL;
//     double t,p,x,y;
//    char *parms[] = { "proj=lcca", "ellps=WGS84", "lon_0=5W", "lat_0=43°N"};
/* *----------------------------------------------------------------------------
    Lambert III IGN*/
//    char *parms[] = { "proj=lcc","ellps=WGS84", "lon_0=2.33722916666666666666E", "lat_0=42.30N", "lat_1=43.2N", "lat_2=45°N"};
/* *----------------------------------------------------------------------------
    Lambert I IGN*/
//    char *parms[] = { "proj=lcc","ellps=WGS84", "lon_0=2.33722916666666666666E", "lat_0=47.70N", "lat_1=°48.59852277777777777777N", "lat_2=50.39591166666666666666N"};
//    char *parms[] = { "proj=merc", "ellps=WGS84", "lon_0=0E", "lat_0=49N", "lat_ts=49N" };
//     const char *parms = "+proj=merc +ellps=WGS84 +lon_0=0E +lat_0=49N +lat_ts=49N";

/* *--------------------------------------------------------------------------------------------------------------------------------------------------
 Projection  		Système  	Parallèle 1  	Parallèle 2  	Latitude origine  	Longitude origine  	X0  		Y0
 Lambert I 		NTF 	48°35'54.682" N 	50°23'45.282" N 	55 gr N 	2°20'14.025" E 		600 000 	1 200 000
 Lambert II 		NTF 	45°53'56.108" N 	47°41'45.652" N 	52 gr N 	2°20'14.025" E 		600 000 	2 200 000
 Lambert II étendu 	NTF 	45°53'56.108" N 	47°41'45.652" N 	52 gr N 	2°20'14.025" E 		600 000 	2 200 000
 Lambert III 		NTF 	43°11'57.449" N 	44°59'45.938" N 	49 gr N 	2°20'14.025" E 		600 000 	3 200 000
 Lambert IV 		NTF 	41°33'37.396" N 	42°46'03.588" N 	43.85 gr N 	2°20'14.025" E 		234.358 	4 185 861.369
 Lambert 93 		RGF93 	44°00'00.000" N 	49°00'00.000" N 	46°30' N 	3°00'00.000" E 		700 000 	6 600 000

+proj=lcc   +lat_1=Latitude of first standard parallel
              +lat_2=Latitude of second standard parallel
              +lat_0=Latitude of false origin
              +lon_0=Longitude of false origin
              +x_0=False Origin Easting
              +y_0=False Origin Northing

proj=lcc  lat_1=  lat_2= lat_0= lon_0= x_0= y_0=

Lambert 93 proj="proj=lcc lat_1=44  lat_2=49 lat_0=46.5 lon_0=3 x_0=700000 y_0=6600000 ellps=WGS84"

Lambert II proj="proj=lcc +lon_0=2.337229166667 +lat_0=46.8 +lat_1=45.898918964419 +lat_2=47.696014502038 +x_0=600000 +y_0=2200000 +a=6378249.2 +rf=293.4660213 +towgs84=-168,-60,320,0,0,0,0 +units=m +no_defs"

Lambert III +proj="proj=lcc +ellps=clrk80 +x_0=600000 +y_0=3200000 +lon_0=2.337291667 +lat_0=44.1 +lat_1=43.19929139 +lat_2=44.99609361"
            +proj="proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=600000 +y_0=3200000 +lon_0=2d20'14.025 +lat_0=44d06' +lat_1=43d11'57.44859 +lat_2=44d59'45.93773"
            
            
proj1=LAMB1 "Lambert I Nord" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=600000 +y_0=200000 +lon_0=2d20'14.025 +lat_0=49d30' +lat_1=48d35'54.682 +lat_2=50d23'45.282"
proj2=LAMB2 "Lambert II Centre" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=600000 +y_0=200000 +lon_0=2d20'14.025 +lat_0=46d48' +lat_1=45d53'56.108 +lat_2=47d41'45.652"
proj3=LAMB3 "Lambert III Sud" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=600000 +y_0=200000 +lon_0=2d20'14.025 +lat_0=44d06' +lat_1=43d11'57.44859 +lat_2=44d59'45.93773"
proj4=LAMB4 "Lambert IV Corse" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=234.358 +y_0=185861.369 +lon_0=2d20'14.025 +lat_0=42d09'54 +lat_1=41d33'37.396 +lat_2=42d46'03.588"
proj5=LAMB93 "Lambert 93" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=700000 +y_0=6600000 +lon_0=3d +lat_0=46d30' +lat_1=44d +lat_2=49d"
proj6=LAMB1C "Lambert I Carto" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=600000 +y_0=1200000 +lon_0=2d20'14.025 +lat_0=49d30' +lat_1=48d35'54.682 +lat_2=50d23'45.282"
proj7=LAMB2C "Lambert II Carto" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=600000 +y_0=2200000 +lon_0=2d20'14.025 +lat_0=46d48' +lat_1=45d53'56.108 +lat_2=47d41'45.652"
proj8=LAMB3C "Lambert III Carto" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=600000 +y_0=3200000 +lon_0=2d20'14.025 +lat_0=44d06' +lat_1=43d11'57.44859 +lat_2=44d59'45.93773"
proj3=LAMB3    "Lambert III Sud" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=600000 +y_0=200000  +lon_0=2d20'14.025 +lat_0=44d06' +lat_1=43d11'57.44859 +lat_2=44d59'45.93773"
proj9=LAMB4C "Lambert IV Carto" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=234.358 +y_0=4185861.369 +lon_0=2d20'14.025 +lat_0=42d09'54 +lat_1=41d33'37.396 +lat_2=42d46'03.588"
proj10=LAMBE "Lambert II Etendu" "+proj=lcc +ellps=clrk80IGN +towgs84=0,0,0 +x_0=600000 +y_0=2200000 +lon_0=2d20'14.025 +lat_0=46d48' +lat_1=45d53'56.108 +lat_2=47d41'45.652"
proj11=RGF93CC42 "Lambert 93 Zone 1" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=1700000 +y_0=1200000 +lon_0=3d +lat_0=42d +lat_1=41d15' +lat_2=42d45'"
proj12=RGF93CC43 "Lambert 93 Zone 2" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=1700000 +y_0=2200000 +lon_0=3d +lat_0=43d +lat_1=42d15' +lat_2=43d45'"
proj13=RGF93CC44 "Lambert 93 Zone 3" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=1700000 +y_0=3200000 +lon_0=3d +lat_0=44d +lat_1=43d15' +lat_2=44d45'"
proj14=RGF93CC45 "Lambert 93 Zone 4" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=1700000 +y_0=4200000 +lon_0=3d +lat_0=45d +lat_1=44d15' +lat_2=45d45'"
proj15=RGF93CC46 "Lambert 93 Zone 5" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=1700000 +y_0=5200000 +lon_0=3d +lat_0=46d +lat_1=45d15' +lat_2=46d45'"
proj16=RGF93CC47 "Lambert 93 Zone 6" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=1700000 +y_0=6200000 +lon_0=3d +lat_0=47d +lat_1=46d15' +lat_2=47d45'"
proj17=RGF93CC48 "Lambert 93 Zone 7" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=1700000 +y_0=7200000 +lon_0=3d +lat_0=48d +lat_1=47d15' +lat_2=48d45'"
proj18=RGF93CC49 "Lambert 93 Zone 8" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=1700000 +y_0=8200000 +lon_0=3d +lat_0=49d +lat_1=48d15' +lat_2=49d45'"
proj19=RGF93CC50 "Lambert 93 Zone 9" "+proj=lcc +ellps=GRS80 +grids=c:\topocad\ign.dat +x_0=1700000 +y_0=9200000 +lon_0=3d +lat_0=50d +lat_1=49d15' +lat_2=50d45'"

27563, "epsg", 27563, "NTF (Paris) / Lambert Sud France",
-        "+proj=lcc +lat_1=44.10000000000001 +lat_0=44.10000000000001 +lon_0=2d20 +k_0=0.999877499 +x_0=600000 +y_0=200000 +a=6378249.2 +b=6356515 +towgs84=-168,-60,320,0,0,0,0 +pm=paris +units=m +no_defs"},

----------------------------------------------------------------------------------------------------------------------------------------------------*/

    /*
    if ( ! (ref = pj_init_plus(parms)) ) TRAP_ERR_EXIT(1, "Projection initialization failed\n");
    pj_free(ref);
    */
    
    double *xx=new double[mesh->nvtxs];
    double *yy=new double[mesh->nvtxs];
    double *tt=new double[mesh->nvtxs];
    double *pp=new double[mesh->nvtxs];
    for(i=0;i<mesh->nvtxs;i++) {
//       xx[i]=mesh->vertices[i].lon+800000;
//       yy[i]=mesh->vertices[i].lat+1800000;
      xx[i]=mesh->vertices[i].lon;
      yy[i]=mesh->vertices[i].lat;
      }
    
    if(proj4!=0) {
      status=projection_to_geo (proj4, xx, yy, mesh->nvtxs);
      for(i=0;i<mesh->nvtxs;i++) {
        mesh->vertices[i].lon=xx[i];
        mesh->vertices[i].lat=yy[i];
        }
      }
    else {
      CONV_LAMBERT_DEGDEC_TAB(mesh->nvtxs,xx, yy, tt, pp );
      for(i=0;i<mesh->nvtxs;i++) {
        mesh->vertices[i].lon=tt[i];
        mesh->vertices[i].lat=pp[i];
        }
      }
    
    delete[] xx;
    delete[] yy;
    delete[] tt;
    delete[] pp;
//     for(i=0;i<mesh->nvtxs;i++) {
//       x=mesh->vertices[i].lon;
//       y=mesh->vertices[i].lat;
//       projection_to_geo(ref,&p,&t,x,y);
//       mesh->vertices[i].lon=t;
//       mesh->vertices[i].lat=p;
//       }
    }
    
  status= fe_edgetable(mesh,0,0);
//   status= discretisation_init(mesh, LGP1,0);
//   status= archiving_UGdummy2D("test.nc", *mesh, "z", "m", z, 0.0, 0, LGP1);
//   status= archiving_UGdummy2D("test.nc", *mesh, "u", "m/s", u, 0.0, 0, LGP1);
//   status= archiving_UGdummy2D("test.nc", *mesh, "v", "m/s", v, 0.0, 0, LGP1);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_read_TELEMAC_BINARY (const char *filename, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  int i,j,status;
  char *title=NULL,*name=NULL,**varname=NULL;
  int recsize,nval;
  int nframes;
  int *ibuf=NULL;
  float *rbuf=NULL;
  int *nbv=NULL;
  bool need_swap=0;
  double *time;
  long pos, first;
  date_t start;
  int nitems;
  
  in = fopen(filename, "r");
  if(in==0) {
    return(-1);
    }

/* *----------------------------------------------------------------------------
  Title */
  fread(&recsize,sizeof(int),1,in);
  
  if(recsize<0)   need_swap=1;
  if(recsize>256) need_swap=1;
  
  if(need_swap) {
    printf("TELEMAC_BINARY, binary will be swapped\n");
    }

  if(need_swap) {
    recsize=lswap(recsize);
    }

  title=new char[recsize+1];
  fread(title,recsize,1,in);
  title[recsize]=0;
  
  printf("TELEMAC_BINARY, title: %s (%d)\n",title,recsize);

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

/* *----------------------------------------------------------------------------
  NBV */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(int);
  nbv=new int[nval];
  printf("TELEMAC_BINARY, variables nrecords : %d\n",nval);

  fread(nbv,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) nbv[i]=lswap(nbv[i]);
    }
  for(i=0;i<nval;i++) printf("TELEMAC_BINARY, variable %d  : %d\n",i,nbv[i]);


  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

/* *----------------------------------------------------------------------------
  NAMES */
  varname=new char*[nbv[0]];
  for(i=0;i<nbv[0];i++) {
    fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }

    name=new char[recsize+1];
    fread(name,recsize,1,in);
    name[recsize]=0;
    
    printf("TELEMAC_BINARY, variable name %d: %s (%d)\n",i,name,recsize);
    varname[i]=strdup(name);
    
    delete[] name;

    fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }
    }

/* *----------------------------------------------------------------------------
  IPARAM */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(int);
  ibuf=new int[nval];
  printf("TELEMAC_BINARY, parameters nrecords : %d\n",nval);

  fread(ibuf,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
    }

  for(i=0;i<nval;i++) printf("TELEMAC_BINARY, parameter %d  : %d\n", i, ibuf[i]);

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  int date=(ibuf[9]==1);

  delete[] ibuf;
  
  if(date==1) {
    fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }

    nval=recsize/sizeof(int);
    ibuf=new int[nval];
    printf("TELEMAC_BINARY, date nrecords: %d\n",nval);

    fread(ibuf,recsize,1,in);
    if(need_swap) {
      for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
      }

    for(i=0;i<nval;i++) printf("TELEMAC_BINARY, date record %d  : %d\n",i,ibuf[i]);

    fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }

    start.year  =ibuf[0];
    start.month =ibuf[1];
    start.day   =ibuf[2];
    start.second=ibuf[3]*3600+ibuf[4]*60+ibuf[5];
  
    delete[] ibuf;
    }
  
/* *----------------------------------------------------------------------------
  CARDINAL */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  nval=recsize/sizeof(int);
  ibuf=new int[nval];
  printf("TELEMAC_BINARY, cardinal nrecords : %d\n",nval);

  fread(ibuf,recsize,1,in);
  if(need_swap) {
    for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
    }
  for(i=0;i<nval;i++) printf("TELEMAC_BINARY, cardinal %d  : %d\n",i,ibuf[i]);

  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  delete[] ibuf;

/* *----------------------------------------------------------------------------
  CONNECTIVITY */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }
  pos=recsize+4;
  fseek(in,pos,SEEK_CUR);
//   nval=recsize/sizeof(int);
//   ibuf=new int[nval];
// 
//   fread(ibuf,recsize,1,in);
//   if(need_swap) {
//     for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
//     }
// 
//   fread(&recsize,sizeof(int),1,in);
//   if(need_swap) {
//     recsize=lswap(recsize);
//     }
// 
//   delete[] ibuf;

/* *----------------------------------------------------------------------------
  BOUNDARY CODES */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  pos=recsize+4;
  fseek(in,pos,SEEK_CUR);
//   nval=recsize/sizeof(int);
//   ibuf=new int[nval];
// 
//   fread(ibuf,recsize,1,in);
//   if(need_swap) {
//     for(i=0;i<nval;i++) ibuf[i]=lswap(ibuf[i]);
//     }
// 
//   fread(&recsize,sizeof(int),1,in);
//   if(need_swap) {
//     recsize=lswap(recsize);
//     }
// 
//   delete[] ibuf;

/* *----------------------------------------------------------------------------
  LONGITUDES */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  pos=recsize+4;
  fseek(in,pos,SEEK_CUR);
//   nval=recsize/sizeof(float);
//   rbuf=new float[nval];
// 
//   fread(rbuf,recsize,1,in);
//   if(need_swap) {
//     for(i=0;i<nval;i++) rbuf[i]=fswap(rbuf[i]);
//     }
// 
//   fread(&recsize,sizeof(int),1,in);
//   if(need_swap) {
//     recsize=lswap(recsize);
//     }
//   delete[] rbuf;

/* *----------------------------------------------------------------------------
  LATITUDES */
  fread(&recsize,sizeof(int),1,in);
  if(need_swap) {
    recsize=lswap(recsize);
    }

  pos=recsize+4;
  fseek(in,pos,SEEK_CUR);
//   nval=recsize/sizeof(float);
//   rbuf=new float[nval];
// 
//   fread(rbuf,recsize,1,in);
//   if(need_swap) {
//     for(i=0;i<nval;i++) rbuf[i]=fswap(rbuf[i]);
//     }
// 
//   fread(&recsize,sizeof(int),1,in);
//   if(need_swap) {
//     recsize=lswap(recsize);
//     }
//   delete[] rbuf;

  float *h=new float[mesh.nvtxs];
  float *u=new float[mesh.nvtxs];
  float *v=new float[mesh.nvtxs];
  float *z=new float[mesh.nvtxs];
  
  for(int n=0;n<mesh.nvtxs;n++) {
    h[n]=0;
    z[n]=0;
    u[n]=0;
    v[n]=0;
    }
  
  FILE *sample=fopen("sample.dat","w");
  float mask=-999.9;
  
  first=ftell(in);
  
  nframes=0;
  while(true) {
/* *----------------------------------------------------------------------------
    TIME */
    nitems=fread(&recsize,sizeof(int),1,in);
    if(nitems!=1) break;
  
    if(need_swap) {
      recsize=lswap(recsize);
      }
    nval=recsize/sizeof(float);
    if(nval==0) {
      printf("empty time record at frame %d\n",nframes);
//      break;
      recsize=4;
      nval=1;
      }
    pos=recsize+4+nbv[0]*(mesh.nvtxs+2)*4;
    status=fseek(in,pos,SEEK_CUR);
    if(status==-1) break;
    nframes++;
    }
    
  time=new double[nframes];
  double tprev;
  
  status= discretisation_init(&mesh, LGP1,0);

  fseek(in,first,SEEK_SET);
  
  nframes=0;
  while(true) {
/* *----------------------------------------------------------------------------
    TIME */
    nitems=fread(&recsize,sizeof(int),1,in);
    if(nitems!=1) break;
  
    if(need_swap) {
      recsize=lswap(recsize);
      }

    nval=recsize/sizeof(float);
    printf("TELEMAC_BINARY, time nrecords: %d %d\n",nval,nframes);
    if(nval==0) {
      printf("empty time record at frame %d\n",nframes);
//      break;
      recsize=4;
      nval=1;
      }
    rbuf=new float[nval];

    fread(rbuf,recsize,1,in);
    if(need_swap) {
      for(i=0;i<nval;i++) rbuf[i]=fswap(rbuf[i]);
      }
    
    double t0=cnes_time(start, 'd');
  
    for(i=0;i<nval;i++) {
      printf("TELEMAC_BINARY, time %d  : %lf\n",nframes,rbuf[i]);
      time[nframes]=t0+rbuf[i]/(24.*3600.);
      if((rbuf[i]==0.)&&(nframes>0)) time[nframes]=tprev+1800/(24.*3600.);
      tprev=time[nframes];
      }
    
    fread(&recsize,sizeof(int),1,in);
    if(need_swap) {
      recsize=lswap(recsize);
      }
    delete[] rbuf;
  
/* *----------------------------------------------------------------------------
    VARIABLES */
    for(j=0;j<nbv[0];j++) {
      fread(&recsize,sizeof(int),1,in);
      if(need_swap) {
        recsize=lswap(recsize);
        }

      nval=recsize/sizeof(float);
      printf("TELEMAC_BINARY, variable %d nrecords: %d\n",j,nval);
      if(nval==0) {
        long pos;
        pos=ftell(in);
        printf("empty record at variable %d, frame %d (pos=%ld bytes)\n",j,nframes, pos);
//        break;
        recsize=4*mesh.nvtxs;
        nval=mesh.nvtxs;
        }

      rbuf=new float[nval];
      
      fread(rbuf,recsize,1,in);
      if(need_swap) {
        for(i=0;i<nval;i++) rbuf[i]=fswap(rbuf[i]);
        }

      fread(&recsize,sizeof(int),1,in);
      if(need_swap) {
        recsize=lswap(recsize);
        }

      if(strstr(varname[j],"SURFACE LIBRE")!=0) {
        for(i=0;i<mesh.nvtxs;i++) {
          z[i]=rbuf[i];
          }
        }
        
      if(strstr(varname[j],"FOND")!=0 or strstr(varname[j],"BOTTOM")!=0) {
        for(i=0;i<mesh.nvtxs;i++) {
          h[i]=rbuf[i];
          }
        }
        
      if(strstr(varname[j],"VITESSE U")!=0) {
        for(i=0;i<mesh.nvtxs;i++) {
          u[i]=rbuf[i];
          }
        }
      if(strstr(varname[j],"VITESSE V")!=0) {
        for(i=0;i<mesh.nvtxs;i++) {
          v[i]=rbuf[i];
          }
        }
      delete[] rbuf;
      }
    i=0;
    fprintf(sample,"%12.6lf %9.4f %9.4f %9.4f\n",time[nframes],z[i],u[i],v[i]);

    status= archiving_UGdummy2D("test.nc", mesh, "z", "m",   z, mask, nframes, LGP1);
    status= archiving_UGdummy2D("test.nc", mesh, "u", "m/s", u, mask, nframes, LGP1);
    status= archiving_UGdummy2D("test.nc", mesh, "v", "m/s", v, mask, nframes, LGP1);
    status= archiving_UGdummy2D("test.nc", mesh, "h", "m/s", h, mask, nframes, LGP1);
    
    status=poc_writetime("test.nc", nframes, "time", time[nframes]*24.*3600.);


    nframes++;
  
    }

  fclose(sample);
  
  fclose(in);
  
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readmesh_MODULEFP2 (const char *filename, mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  int i,m,nnodes,nn,np;
  int nitems,n1,n2,n3,n4,n5,n6,count,status;
  char *s=NULL;
  char line[1024],separator[32];
  numerics_t LGP2_C_numerics;

/*------------------------------------------------------------------------------
  element and node list format*/
  in = fopen(filename, "r");
  if(in==0)
      return(-1);

/*------------------------------------------------------------------------------
  rough estimate of the maximum count of nodes*/
  count=0;
  do {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"(NTRI ) :");
    } while(s==0);
  sscanf((s+9),"%d", &(mesh->ntriangles));
  mesh->triangles= new triangle_t [mesh->ntriangles];
  if(mesh->triangles == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }

  do {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"(NOE  ) :");
    } while(s==0);
  sscanf((s+9),"%d", &nnodes);

  LGP2_C_numerics.nnodes=nnodes;
  LGP2_C_numerics.type=LGP2_C;
  LGP2_C_numerics.nodes=new int*[mesh->ntriangles];
  for(m=0;m<mesh->ntriangles;m++) {
    LGP2_C_numerics.nodes[m]=new int[6];
    }

  do {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"(NP  ) :");
    } while(s==0);
  sscanf((s+9),"%d", &(mesh->nvtxs));
  mesh->vertices= new vertex_t [mesh->nvtxs];
  if(mesh->vertices == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }

  do {
    fgets(line,1024,in);
    if(feof(in)) break;
    s=strstr(line,"| POINT |");
    } while(s==0);
  fgets(line,1024,in);

  for(i=0;i<mesh->nvtxs;i++) {
    nitems=fscanf(in, "%s %d %s %lf %s %lf %s", &separator, &n1,
                                                &separator, &mesh->vertices[i].lon,
                                                &separator, &mesh->vertices[i].lat,
                                                &separator);
    }

  for(m=0;m<mesh->ntriangles;m++) {
    do {
      fgets(line,1024,in);
      if(feof(in)) break;
      s=strstr(line,"NOMBRE DE  NOEUDS :");
      } while(s==0);
//    sscanf((s+19),"%d", &nn);
    sscanf((s+19),"%d %d %d %d %d %d %d", &nn,&n1,&n3,&n5,&n2,&n4,&n6);
    LGP2_C_numerics.nodes[m][0]=n1-1;
    LGP2_C_numerics.nodes[m][1]=n2-1;
    LGP2_C_numerics.nodes[m][2]=n3-1;
    LGP2_C_numerics.nodes[m][3]=n4-1;
    LGP2_C_numerics.nodes[m][4]=n5-1;
    LGP2_C_numerics.nodes[m][5]=n6-1;
    do {
      fgets(line,1024,in);
      if(feof(in)) break;
      s=strstr(line,"NOMBRE DE  POINTS :");
      } while(s==0);
    sscanf((s+19),"%d %d %d %d", &np,&n1,&n2,&n3);
    mesh->triangles[m].vertex[0]=n1-1;
    mesh->triangles[m].vertex[1]=n2-1;
    mesh->triangles[m].vertex[2]=n3-1;
    }

  status=fe_e2n(mesh);

  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int parse_MeshOptions(const string options, char* & ElementFile, char* & NodeFile, char* & DepthFile)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<string> tokens,keys,values;
  int k;
  string delimiter=" ";
  
  ElementFile=NodeFile=DepthFile=0;
  
  delimiter=" ";
  tokens=string_split(options, delimiter);
  
  delimiter="=";
  for(k=0;k<tokens.size();k++) {
    vector<string> tmp=string_split(tokens[k], delimiter);
    keys.push_back(tmp[0]);
    values.push_back(tmp[1]);
    }
    
  for(k=0;k<keys.size();k++) {
    if(keys[k]=="elements") {
      ElementFile=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="nodes") {
      NodeFile=strdup(values[k].c_str());
      continue;
      }
    if(keys[k]=="depths") {
      DepthFile=strdup(values[k].c_str());
      continue;
      }
    }
  return(0);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readmesh(const char *filename, int fmt, mesh_t *mesh, const char *proj4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  char *element_filename=NULL, *node_filename=NULL, *depth_filename=NULL;
  char *rootname=NULL;
  extern int fe_read_TELEMAC_BINARY (const char *filename, mesh_t & mesh);
  bool debug=true;
  
  if(fmt==MESH_FILE_FORMAT_UNKNOWN)
    fmt=fe_find_format(filename);
  
  switch(fmt) {
    case MESH_FILE_FORMAT_TRIGRID:
/*------------------------------------------------------------------------------
      neighbour list format*/
      status=fe_readmesh_NGH(filename, mesh);
      break;

    case MESH_FILE_FORMAT_TRIANGLE:
      status=fe_readmesh_TGL(element_filename, node_filename, mesh);
      break;

    case MESH_FILE_FORMAT_QUODDY:
      status=parse_MeshOptions((string) filename, element_filename, node_filename, depth_filename);
      status=fe_readmesh_QUODDY(element_filename, node_filename, depth_filename, mesh, debug);
      break;

    case MESH_FILE_FORMAT_GOM:
      status=-1;
      break;

    case MESH_FILE_FORMAT_NC3D:
      status=fe_readmesh3d(filename, mesh, 0);
      break;

    case MESH_FILE_FORMAT_MODULEF_P1:
      status=-1;
      break;

    case MESH_FILE_FORMAT_MODULEF_P2:
      status=fe_readmesh_MODULEFP2 (filename, mesh);
      break;

    case MESH_FILE_FORMAT_GMSH:
      status=fe_readmesh_GMSH(filename, mesh);
      break;
      
    case MESH_FILE_FORMAT_GMSH_WW:
      status=fe_readmesh_GMSH_WW(filename, mesh);
      break;

    case MESH_FILE_FORMAT_SCHISM:
      status=fe_readmesh_SCHISM(filename, mesh);
      break;

    case MESH_FILE_FORMAT_TELEMAC_ASCII:
      rootname=strdup(filename);
      rootname[strlen(filename)-4]=0;
      status=fe_readmesh_TELEMAC(rootname, mesh);
      break;

    case MESH_FILE_FORMAT_TELEMAC_BINARY:
      status=fe_readmesh_TELEMAC_BINARY(filename, mesh, proj4);
      status=fe_read_TELEMAC_BINARY (filename, *mesh);
      break;

    default:
      status=-1;
      break;

    }
    
  if(status) TRAP_ERR_RETURN(status,1,"could not load mesh from %s (format=%d)\n", filename,fmt);
  
  status=fe_initaffine(mesh);
  
  return(status);
}


/*----------------------------------------------------------------------------*/
/** reads the boundary code from a belfile and save it into a mesh_t object (mesh.edges.code)
\date reviewed 1 Aug 2011
\author Damien Allain

\todo more documentation

@param belfile the path+filename of the belfile to be read
@param mesh
\param level How to behave when the belfile can not be opened. Default:0.
  0 : do not give any message, 1 and 2: give a message
  0 and 1 : returns -1, 2: exits with an error code
\returns 0 on success, -1 if error
*/
/*----------------------------------------------------------------------------*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void read_boundarycode_01(const char *filename, mesh_t & mesh, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------

  construct boundary codes from historical BEL file

  Input:
   -  mesh.vertices cross tables needed
   -  mesh.triangles cross tables needed

  assigns the boundary codes as follow (see also mesh.def):
 
  MESH_UNDEFINED_NODE     -1
  MESH_INTERIOR_NODE       0
  MESH_LAND_NODE           1
  MESH_ISLAND_NODE         2
  MESH_PERMEABLE_NODE      3
  MESH_GEOSTROPHY_NODE     4
  MESH_ELEVATION_NODE      5
  MESH_PERMEABLE_NODE      6   rivers with elevation (amazon...)

----------------------------------------------------------------------*/
{
  int dum1, dum2, n, nd1, nd2, code;
  int nitems, found;
  int nedges, count = 0;
  char char_dum[500];
  FILE *in;
  edge_t *edges;

  nedges = 0;

  edges  = mesh.edges;
  nedges = mesh.nedges;

/**----------------------------------------------------------------------
  initialize edges' code to UNDEFINED VALUE */
  for(n = 0; n < nedges; n++)
    edges[n].code = MESH_UNDEFINED_NODE;

  in = fopen(filename, "r");
  if(in==0) {
    TRAP_ERR_EXIT(-1,"read_boundarycode failed (cannot open file: %s)",filename);
    }

  fgets(char_dum, n, in);
  fgets(char_dum, n, in);

  count = 0;

  while(!feof(in)) {
    nitems = fscanf(in, "%d %d %d %d %d", &dum1, &nd1, &nd2, &dum2, &code);
    if(nitems != 5)
      goto done;
/**----------------------------------------------------------------------
    C-array from fortan numbering*/
    nd1--;
    nd2--;
    count++;
/**----------------------------------------------------------------------
    find corresponding triangle edge */
    found = 0;
    n=fe_isedge(mesh, nd1,nd2);
    if(n!=-1) {
      edges[n].code=code;
      found=1;
      }
    if(found != 1) {
      TRAP_ERR_EXIT(-1,"boundary error at edge %d, nodes %d %d, code %d --- edge not found\n", dum1, nd1, nd2, code);
      }
    }

done:

  fclose(in);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void write_boundarycode_01(const char *filename, mesh_t & mesh, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------

  Write boundary codes in the historical BEL file format
  
----------------------------------------------------------------------*/
{
  int k, l;
  int dum, n, nd1, nd2;
  int nitems;
  int count = 1;
  FILE *in;

/**----------------------------------------------------------------------
   */
  in = fopen(filename, "w");
  if(in==0) {
    TRAP_ERR_EXIT(-1,"read_boundarycode failed (cannot open file: %s)",filename);
    }

  fprintf(in,"automatically generated by tools softwares\n");
  fprintf(in,"automatically generated by tools softwares\n");
  
  dum=0;

  for(l = 0; l < mesh.nlimits; l++) {
    for(k = 0; k < mesh.limits[l].nedges; k++) {
      n=mesh.limits[l].edges[k];
      nd1=mesh.edges[n].extremity[0];
      nd2=mesh.edges[n].extremity[1];
      nitems = fprintf(in, "%6d %6d %6d %2d %2d\n", count, nd1+1, nd2+1, dum, mesh.edges[n].code);
      count++;
      }
    }

  count = 0;

  fclose(in);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void read_boundarycode_02(const char *filename, mesh_t & mesh, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------

  construct boundary codes from new BEL file format

  Input:
   -  mesh.vertices cross tables needed
   -  mesh.triangles cross tables needed

  MESH_UNDEFINED_NODE     -1
  MESH_INTERIOR_NODE       0
  MESH_LAND_NODE           1
  MESH_ISLAND_NODE         2
  MESH_PERMEABLE_NODE      3
  MESH_GEOSTROPHY_NODE     4
  MESH_ELEVATION_NODE      5

----------------------------------------------------------------------*/
{
  int k, l, status;
  int n, n1, n2;
  int nitems;
  int code, count = 0;
  char char_dum[1024];
  FILE *in;
//  edge_t *edges;
  double x1,y1,x2,y2;

/**----------------------------------------------------------------------
  initialize edges' code to UNDEFINED VALUE */
  for(n = 0; n < mesh.nedges; n++)
    mesh.edges[n].code = MESH_UNDEFINED_NODE;
  
  for(l = 0; l < mesh.nlimits; l++) {
    for(k = 0; k < mesh.limits[l].nedges; k++) {
      n=mesh.limits[l].edges[k];
      mesh.edges[n].code = 1;
      }
    }

  in = fopen(filename, "r");
  if(in == NULL) TRAP_ERR_EXIT(errno,"%s(\"%s\",,) failed (%s)\n",__func__,filename,strerror(errno));

  n=1024;
  fgets(char_dum, n, in);

  count = 0;

  while(!feof(in)) {
    char *s=fgets(char_dum, n, in);
    if(s==0) break;
    if(s[0]=='#') continue;
    nitems = sscanf(char_dum, "%lf %lf %lf %lf %d", &x1, &y1, &x2, &y2, &code);
    if(code!=5) {
      TRAP_ERR_EXIT(errno,"%s(\"%s\",,) failed (%s)\n",__func__,filename,strerror(errno));
      }
    if(nitems < 4)
      goto done;
    if(nitems == 4)
      code=MESH_ELEVATION_NODE;
    n1=fe_nearest_boundaryvertex (mesh,x1,y1,0,0);
    n2=fe_nearest_boundaryvertex (mesh,x2,y2,0,0);
    status=fe_setOBCflags(mesh, n1, n2, code);
    }

done:

  fclose(in);
  write_boundarycode_01("old-format.bel", mesh, level);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_read_boundarycode(const char *filename, mesh_t & mesh, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  char dum[500];
  int fmt;
  
  in = fopen(filename, "r");
  if(in == NULL) TRAP_ERR_EXIT(errno,"%s(\"%s\",,) failed (%s)\n",__func__,filename,strerror(errno));
  
  fgets(dum, 500, in);
  char *s=strstr( (char *) dum, "#OPEN BOUNDARIES");
  if(s==0) {
    fmt=0;
    }
  else {
    fmt=1;
    }
  fclose(in);
  
  switch (fmt) {
    case 0:
      read_boundarycode_01(filename, mesh, level);
      break;
    case 1:
      read_boundarycode_02(filename, mesh, level);
      break;
  }
  
  return(0);
}


/*----------------------------------------------------------------------------*/
/**
\date reviewed 1 Aug 2011
\author Damien Allain

\todo more documentation (when necessary)

\sa fe_read_boundarycode().

@param belfile the path+filename of the belfile to be read
@param mesh
\param level How to behave when the belfile can not be opened. Default:0.
  0 : do not give any message, 1 and 2: give a message
  0 and 1 : returns -1, 2: exits with an error code
\returns 0 on success, -1 if error
*/
/*----------------------------------------------------------------------------*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_write_boundarycode(const char *belfile, const mesh_t &mesh, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  dum=0,n, n1, n2, k, l;
  int  nitems,count=0;
  FILE *in=NULL;

  in = fopen(belfile, "w");
  if (in == NULL) {
    if(level){
      check_error(errno,__LINE__,__FILE__,"fopen error with %s",belfile);
      }
    if(level>=2) TRAP_ERR_EXIT(-1,"exiting\n");
    return -1;
    }

  fprintf(in,"%s\n","boundary codes");
  fprintf(in,"%s\n","comments");

  count=1;
  for(l = 0; l < mesh.nlimits; l++) {
    for(k = 0; k < mesh.limits[l].nedges; k++) {
      n=mesh.limits[l].edges[k];
      n1= mesh.limits[l].vertex[k]+1;
      n2= mesh.limits[l].vertex[(k+1)%mesh.limits[l].nvertex] + 1;
#if (0)
      nitems=fprintf(in, "%d %d %d %d %d\n", count, n1, n2, dum, 1);
#else
     nitems=fprintf(in, "%d %d %d %d %d\n", count, n1, n2, dum, mesh.edges[n].code);
#endif
      count++;
      }
    }
  fclose(in);

  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_write_boundarycode(const char *belfile, vector<plg_t> & boundaries)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int last, status;
  FILE *in=NULL;
  vector<plg_t> polygons, opened;
  
  in = fopen(belfile, "w");
  if (in == NULL) STDERR_BASE_LINE("fopen error with %s",belfile);
  fprintf(in,"%s\n","#OPEN BOUNDARIES");

  status=plg_setdirect(boundaries);
  
  polygons=plg_split(boundaries,false);
  for(int s=0;s<polygons.size();s++) {
    if(polygons[s].flag[0]=='M' or polygons[s].flag[0]=='X') {
      opened.push_back(polygons[s]);
      }
    }
    
  printf("#################################################################\n");
  printf("create T-UGOm BEL file: %d segments found (%d open)\n", polygons.size(), opened.size());
      
  for(int s=0;s<opened.size();s++) {
    if(plg_isclosed(opened[s])==1) last=opened[s].npt-2;
    else last=opened[s].npt-1;
    fprintf(in,"%lf %lf %lf %lf    5\n",opened[s].t[0], opened[s].p[0],opened[s].t[last], opened[s].p[last]);
    }
  
  fclose(in);
  
  opened.clear();
  polygons.clear();
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_SaveLimits(const char *filename, const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  int k,l,n;

  in = fopen(filename, "w");
  if(in==0)
      return(-1);

  for (l=0; l<mesh.nlimits; l++) {
    fprintf(in, "%d %d \n", l+1, mesh.limits[l].nvertex+1);
    for (k=0; k<mesh.limits[l].nvertex;k++) {
      n=mesh.limits[l].vertex[k];
      fprintf(in, "%d %lf %lf\n", k+1, mesh.vertices[n].lon, mesh.vertices[n].lat);
      }
    n=mesh.limits[l].vertex[0];
    fprintf(in, "%d %lf %lf\n", k+1, mesh.vertices[n].lon, mesh.vertices[n].lat);
    }

  fclose(in);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readnodes_TGD(const char *filename,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------

  TRIGRID format, containing some structure informations:
  - boundary nodes and their sequence
  - interior nodes

------------------------------------------------------------------------*/
{
  FILE *in=NULL;
  vertex_t *set=NULL;
  int i,j,k,nndes,l,ninteriors;

  in = fopen(filename, "r");
  if(in==0)
      return(-1);

  fscanf(in, "%d ", &mesh->nvtxs);
//  printf("#Number of nodes: %d \n",mesh->nvtxs);
//  fscanf(in, "%d ", &nexteriors);
  fscanf(in, "%d ", &(mesh->nlimits));

  nndes=mesh->nvtxs;
  exitIfNull(
    mesh->vertices=new vertex_t[nndes]
    );
  set= mesh->vertices;

  if(set == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }

  for (i=0; i<nndes; i++) {
    set[i].ngh= NULL;
    set[i].nngh= 0;
    }

  i=0;

  if(mesh->nlimits!=0) {
    mesh->limits=new limit_t[mesh->nlimits];
    for (l=0; l<mesh->nlimits; l++) {
      fscanf(in, "%d", &(mesh->limits[l].nvertex));
      mesh->limits[l].vertex=new int[mesh->limits[l].nvertex];
      for (k=0; k<mesh->limits[l].nvertex; k++) {
        fscanf(in,"%lf %lf %f\n",&(set[i].lon),&(set[i].lat),&(set[i].h));
        mesh->limits[l].vertex[k]=i;
        i++;
        }
      }
    }
  fscanf(in, "%d", &ninteriors);
  for (j=0; j<ninteriors; j++) {
    fscanf(in,"%lf %lf %f",&(set[i].lon),&(set[i].lat),&(set[i].h));
//    fscanf(in,"%d %lf %lf %f",&dum,&(set[i].lon),&(set[i].lat),&(set[i].h));
    i++;
    }

  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readnodes_XYZ(const char *filename,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------

  XYZ format, no  structure informations

------------------------------------------------------------------------*/
{
  FILE *in=NULL;
  vertex_t *set=NULL;
  int i,ninteriors,nexteriors;
  int downward=0,cartesian=0,mode=0,complete=0,masked=0;
  char line[1024];

  in = fopen(filename, "r");
  if(in==0)
      return(-1);
  fgets(line, 1024, in);

/* *-----------------------------------------------------------------------------
  decode header informations*/
  mode     =(strstr(line, "YXZ")!=0);
  complete =(strstr(line, "COMPLETE")!=0);
  cartesian=(strstr(line, "CARTESIAN")!=0);
  downward =(strstr(line, "DOWNWARD")!=0);
  masked   =(strstr(line, "MASK=")!=0);

  fscanf(in, "%d ", &mesh->nvtxs);
//  printf("#Number of nodes: %d \n",mesh->nvtxs);

  nexteriors=0;
  ninteriors=mesh->nvtxs;

  mesh->vertices= new vertex_t[mesh->nvtxs];

  set= mesh->vertices;

  if(set == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }

  for (i=0; i<mesh->nvtxs; i++) {
    set[i].ngh= NULL;
    set[i].nngh= 0;
    fscanf(in, "%lf %lf %lf",&set[i].lon,&set[i].lat,&set[i].h);
    }
  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readnodes_TGL(const char *filename,mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  int dummy,c;
  vertex_t *set=NULL;
  int i,nndes;

  in = fopen(filename, "r");
  if(in==0) {
    printf("read error");
    return(-1);
    }

/**----------------------------------------------------------------------------
  on compte le nombre de ligne */
  nndes=0;
  while (!feof(in)) {
    do { c=fgetc(in);   }  while ((c != '\n') && !feof(in));
    nndes++;
    }

  rewind(in);
  nndes=nndes-1;
  mesh->nvtxs=nndes;

  exitIfNull(mesh->vertices=new vertex_t[nndes]);
  set= mesh->vertices;

  if(set == NULL) {
    printf("allocation error vertex set\n");
    return(-1);
    }
  for (i=0; i<nndes; i++) {
  if( fscanf(in,"%d %lf %lf",&dummy,&set[i].lon,&set[i].lat)!=3)
    printf("scanf error\n");
    }

  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readnodes(const char *filename, int fmt, mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
/*------------------------------------------------------------------------------
  neighbour list format*/

  switch(fmt) {
    case NODE_FILE_FORMAT_TRIGRID:
      status=fe_readnodes_TGD(filename,mesh);
      break;

    case NODE_FILE_FORMAT_TRIANGLE:
      status=-1;
      break;

    case NODE_FILE_FORMAT_QUODDY:
      status=-1;
      break;

//    case NODE_FILE_FORMAT_GOM:
      status=-1;
      break;

//    case NODE_FILE_FORMAT_NC3D:
      status=-1;
      break;

//    case NODE_FILE_FORMAT_MODULEF_P1:
      status=-1;
      break;

//    case NODE_FILE_FORMAT_MODULEF_P2:
      status=-1;
      break;

    case NODE_FILE_FORMAT_GMSH:
      status=-1;
      break;

    default:
      status=-1;
      break;

    }
  return(status);
}
