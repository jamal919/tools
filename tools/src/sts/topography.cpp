
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

***************************************************************************/
#include "config.h"
#include "tugo-prototypes.h"
#include "grd.h"
#include "map.h"
#include "datastream.h"

extern  int fe_ascii_loadr1(const char *, mesh_t, float *, float *, int);
extern  int fe_ascii_loadr2(const char *, mesh_t, float *, float *, float *, int);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void topo_gradientLGP1(mesh_t & mesh, parameter_t P1_data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, m, n, status;
  float *h, *N, *r, *dhdx, *dhdy;
  char filename[1024], rootname[1024];
  double depth, gradx, grady;

/*----------------------------------------------------------------------
  bottom topography gradient computation*/
  for(n = 0; n < mesh.nvtxs; n++) {
    depth = -P1_data.h[n];
    gradx = grad_opx[n][0] * depth;
    grady = grad_opy[n][0] * depth;
    for(k = 0; k < mesh.vertices[n].nngh; k++) {
      m = mesh.vertices[n].ngh[k];
      if(m != -1) {
        depth = -P1_data.h[m];
        gradx += grad_opx[n][k] * depth;
        grady += grad_opy[n][k] * depth;
        }
      }
    P1_data.dhdx[n] = gradx;
    P1_data.dhdy[n] = grady;
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void topo_gradientCQP1(mesh_t & mesh, parameter_t P1_data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, m, n, status;
  float *h, *N, *r, *dhdx, *dhdy;
  char filename[1024], rootname[1024];
  double depth, gradx, grady;

/*----------------------------------------------------------------------
  bottom topography gradient computation*/
  for(n = 0; n < mesh.nvtxs; n++) {
    depth = -P1_data.h[n];
    gradx = grad_opx[n][0] * depth;
    grady = grad_opy[n][0] * depth;
    for(k = 0; k < mesh.vertices[n].nngh; k++) {
      m = mesh.vertices[n].ngh[k];
      if(m != -1) {
        depth = -P1_data.h[m];
        gradx += grad_opx[n][k] * depth;
        grady += grad_opy[n][k] * depth;
        }
      }
    P1_data.dhdx[n] = gradx;
    P1_data.dhdy[n] = grady;
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void topo_gradient(mesh_t & mesh, parameter_t VertexData)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
    switch(mesh.nature()) {
    case FE_TRIANGLE:
      topo_gradientLGP1(mesh, VertexData);
      gP1discretisation=LGP1;
      break;
    case FE_QUADRANGLE:
      topo_gradientCQP1(mesh, VertexData);
      gP1discretisation=CQP1;
      break;
    default:
      topo_gradientLGP1(mesh, VertexData);
      gP1discretisation=-1;
      break;
    }

  
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int topo_init(mesh_t & mesh,parameter_t & IO_data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, m, n, n1, n2, status;
  float *zmin,*buffer,mask;
  char filename[1024], rootname[1024];
  double x,y;
  grid_t grid;

#ifdef DRY
  printf("Warning: drying mode activated, minimum depth still applies (%f)\n", gMinDepth);
#endif

/**---------------------------------------------------------------------
  Database read in LGP1 structure (IO_data) and is handled in temporary P1_data*/

  for(n = 0; n < mesh.nedges; n++) {
    if(mesh.edges[n].code==1) {
      n1 = mesh.edges[n].extremity[0];
      n2 = mesh.edges[n].extremity[1];
/*
      if(IO_data.h[n1] > 15) {
//        printf("%d %f \n",n1,P1_data.h[n1]);
        IO_data.h[n1] = 15;
        }
      if(IO_data.h[n2] > 15) {
//        printf("%d %f \n",n2,P1_data.h[n2]);
        IO_data.h[n2] = 15;
        }
*/
      }
    }

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check :

  Notes

  07/01/2010:

    The following was removed:

#ifdef DRY
/**---------------------------------------------------------------------
  Fred: dry addition/change : no mindepth
  return(0);
#endif

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  zmin = new float[mesh.nvtxs];

  if(strcmp(ZminFile, "NONE") != 0) {
/**---------------------------------------------------------------------
    use local minimum elevation to set minimum depth*/
//    strcpy(filename, rootname);
//    strcat(filename, SlopeFile);
    status=grd_loadgrid(ZminFile,&grid);
    if(status != 0)
      check_error(-1, "file input failure", __LINE__, __FILE__, 1);
    buffer = new float[grid.nx*grid.ny];
    grid.modeH=0;
    status=  grd_loadr1(ZminFile,grid,buffer,&mask);
    
    mask=buffer[0];
    for(n = 0; n < mesh.nvtxs; n++) {
      x = IO_data.lon[n] * r2d;//tools update
      y = IO_data.lat[n] * r2d;//tools update
      if(x < grid.xmin - grid.dx / 2.)
        x = x + 360.0;
      if(x > grid.xmax + grid.dx / 2.)
        x = x - 360.0;
      status = map_interpolation(grid, buffer, mask, x, y,&(zmin[n]));
      if(zmin[n]==mask) {
         zmin[n]=gRelativeMinDepth-gMinDepth;
         }
      if(zmin[n]>0) {
         printf("abnormal zmin at vertex %6d: lon=%lf lat=%lf zmin=%g\n",n,x,y,zmin[n]);
         }
      if(IO_data.h[n]+zmin[n] < gRelativeMinDepth) {
        IO_data.h[n] = gRelativeMinDepth-zmin[n];
        }
      }
    delete[] buffer;
    }
  else {

/**---------------------------------------------------------------------
    use overall minimum depth*/
    for(n = 0; n < mesh.nvtxs; n++) {
      if(IO_data.h[n] < gMinDepth) {
        IO_data.h[n] = gMinDepth;
        }
      }
    }

/**---------------------------------------------------------------------
  copy to mesh structure's depth*/
  for(n = 0; n < mesh.nvtxs; n++) {
    mesh.vertices[n].h = IO_data.h[n];
    }

  delete[] zmin;
}

