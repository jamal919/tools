#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.h"
#include "netcdf-proto.h"
#include "zapper.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_setsigma01(mesh_t * mesh, int nlayers, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------
  cosinus distribution*/
  vertex_t *set;
  int i, j, nndes, nlevels, status;
  double *x, range;
  char output[1024];

  nndes = mesh->nvtxs;
  nlevels=nlayers+1;

  mesh->nlayers = nlayers;
  mesh->nlevels = nlevels;

  set = mesh->vertices;

  x = new double[nlevels];

  range = ((double) nlevels - 1.);

  for(j = 0; j < nlevels; j++) {
    x[j] = 0.5*(1.0 + cos(j * M_PI / range));
    }

  for(i = 0; i < nndes; i++) {
    if(set[i].sigma != NULL) zaparr(set[i].sigma);
    set[i].sigma   = new double[nlevels];
    set[i].zlevels  = new double[nlevels];
    for(j = 0; j < nlevels; j++) {
      set[i].sigma[j]  = -x[j];
      set[i].zlevels[j] = -set[i].h * x[j];
      }
    }

  delete[] x;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_setsigma02(mesh_t * mesh, int nlayers, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------
  uniform distribution*/
  vertex_t *set;
  int i, j, nndes, nlevels, status;
  double *x, range;
  char output[1024];

  nndes = mesh->nvtxs;
  nlevels=nlayers+1;

  mesh->nlayers = nlayers;
  mesh->nlevels = nlevels;

  set = mesh->vertices;

  x = new double[nlevels];

  range = ((double) nlevels - 1.);

  for(j = 0; j < nlevels; j++) { /*Warning: TEST !!!*/
    x[j] = 1.0- (double) j / range;
    }

  for(i = 0; i < nndes; i++) {
    if(set[i].sigma != NULL) zaparr(set[i].sigma);
    set[i].sigma   = new double[nlevels];
    set[i].zlevels  = new double[nlevels];
    for(j = 0; j < nlevels; j++) {
      set[i].sigma[j]  = -x[j];
      set[i].zlevels[j] = -set[i].h * x[j];
      }
    }

  delete[] x;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_setsigma03(mesh_t * mesh, int nlayers, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------
  generalized sigma distribution*/
  vertex_t *set;
  int i, j, nndes, nlevels, status;
  double *x,*d, range;
  char output[1024];

  nndes = mesh->nvtxs;
  nlevels=nlayers+1;

  mesh->nlayers = nlayers;
  mesh->nlevels = nlevels;

  set = mesh->vertices;

  d = new double[nlevels];
  x = new double[nlevels];

  range = ((double) nlevels - 1.);

/*------------------------------------------------------------------------
  first step*/
  for(j = 0; j < nlevels; j++) {
    d[j] = ((float) j)*((float) (nlevels-j));
    }

/*------------------------------------------------------------------------
  second step*/
  x[0]=0;
  for(j = 1; j < nlevels; j++) {
    x[j] = x[j-1]+d[j];
    }

/*------------------------------------------------------------------------
  normalization step*/
  x[0]=0;
  for(j = 1; j < nlevels; j++) {
    x[j] = x[j]/x[nlevels-1];
    }

  for(i = 0; i < nndes; i++) {
    if(set[i].sigma != NULL) zaparr(set[i].sigma);
    set[i].sigma   = new double[nlevels];
    set[i].zlevels  = new double[nlevels];
    for(j = 0; j < nlevels; j++) {
      set[i].sigma[j]  =  x[j];
      set[i].zlevels[j] = -set[i].h * x[j];
      }
    }


  delete[] x;
  delete[] d;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_setsigma(mesh_t * mesh, int nlayers, int type)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,mode=-1;
  char output[1024];

  mode=2;

  switch (mode) {
    case 0:
      status=fe_setsigma01(mesh, nlayers, type);
      break;

    case 1:
      status=fe_setsigma02(mesh, nlayers, type);
      break;

    case 2:
      status=fe_setsigma03(mesh, nlayers, type);
      break;
    }

  strcpy(output, "3Dmesh");
  strcat(output, ".nc");

  status = fe_savemeshNC3D(output, *mesh, 1);

 return (0);
}
