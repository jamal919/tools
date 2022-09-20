
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief finite element interpolation functions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>


#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.def"
#include "fe.h"
#include "geo.h"
#include "map.h"
#include "poc-netcdf-data.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_grid_start(const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// gives the grid that overlaps the given mesh with the given size
/*----------------------------------------------------------------------------*/
{
  grid_t grid;
  int n=0;
  double x,y;
  
  grid.xmin=+INFINITY;
  grid.ymin=+INFINITY;
  
  grid.xmax=-INFINITY;
  grid.ymax=-INFINITY;
  
  for (n=0;n<mesh.nvtxs;n++) {
    vertex_t *vertex=&mesh.vertices[n];
    x=vertex->lon;
    y=vertex->lat;
    updatemin(&grid.xmin,x);
    updatemin(&grid.ymin,y);
    updatemax(&grid.xmax,x);
    updatemax(&grid.ymax,y);
    }
  
  grid.nz  = 1;
  
  return(grid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_grid_n(const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// gives the grid that overlaps the given mesh with the given size
/**
\param mesh mesh to overlap
\todo nx and ny parameters (currently only 400)
\sa get_grid_d()
*/
/*----------------------------------------------------------------------------*/
{
  grid_t grid;
  
  grid=get_grid_start(mesh);
  
  grid.nx  = 400;
  grid.ny  = 400;
  grid.modeH  = 0;
  
  grid.dx  = (grid.xmax-grid.xmin)/(grid.nx-1);
  grid.dy  = (grid.ymax-grid.ymin)/(grid.ny-1);
  
  return(grid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_grid_d(const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// gives the grid that overlaps the given mesh at the given resolution
/**
\param mesh mesh to overlap
\todo dx and dy parameters (currently only 0.1)
\sa get_grid_n()
*/
/*----------------------------------------------------------------------------*/
{
  grid_t grid;
  
  grid=get_grid_start(mesh);
  
  grid.dx  = 0.1;
  grid.dy  = 0.1;
  grid.nx  = (int)( (grid.xmax-grid.xmin)/grid.dx+1 );
  grid.ny  = (int)( (grid.ymax-grid.ymin)/grid.dy+1 );
  
  grid.modeH  = 0;
  
  return(grid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_intpl_LGP0_template(const mesh_t & mesh,const T *buffer, T mask, double t, double p, int n, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *z = buffer[n];
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_intpl_LGP0(const mesh_t & mesh, float *buffer, float mask, double t, double p, int n, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_LGP0_template(mesh, buffer, mask, t, p, n, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_intpl_LGP0(const mesh_t & mesh, double *buffer, double mask, double t, double p, int n, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_LGP0_template(mesh, buffer, mask, t, p, n, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_test_area(const mesh_t & mesh,int m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(mesh.triangles[m].Area==0.0) {
    TRAP_ERR_RETURN(MESH_STATUS_FLAT_ELEMENT,1,"%s:triangle %d is flat !\n",__func__,m);
    }
  if(isnan(mesh.triangles[m].Area)) {
    TRAP_ERR_RETURN(MESH_STATUS_FLAT_ELEMENT,1,"%s:triangle %d is WRONG !\n",__func__,m);
    }
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> int fe_intpl_LGP1_template(const mesh_t & mesh,const T *buffer,T mask,double t,double p,int m,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  LGP1 interpolation

  mesh
  *buffer array[mesh.ntriangles]
  mask mask value of buffer
  t longitude
  p latitude
  m triangle index. See fe_detect()
  *z interpolated value. If at least one value of the triangle is masked, will be set to mask
 
 returns 0

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  T zz;
  double x, y, beta[3];
  int status;
  triangle_t & element=mesh.triangles[m];

  status=fe_test_area(mesh,m);
  if(status) {
    TRAP_ERR_RETURN(status,1,"%s:fe_test_area(,%d)=%d\n",__func__,m,status);
    }

  status = fe_affine_directe(element, t, p, &x, &y, mesh.type);
  if(status) TRAP_ERR_RETURN(status,1,"%s:fe_affine_directe() error\n",__func__);

  status = fe_LGPbase(mesh, LGP1, x, y, beta);
  if(status) TRAP_ERR_RETURN(status,1,"%s:fe_LGPbase(,,%g,%g) error on (%g;%g)\n",__func__,x,y,t,p);
  
  zz = 0;
  for(int i = 0; i < 3; i++) {
    const size_t n = element.vertex[i];
    if(buffer[n]==mask){// NAN!=NAN always
      *z=mask;
      return(status);
      }
    zz += buffer[n] * (T) beta[i];
    }

  *z = zz;
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_intpl_LGP1(mesh_t & mesh, float *buffer, float mask, double t, double p, int m, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_intpl_LGP1_template(mesh,buffer,mask,t,p,m,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_intpl_LGP1(mesh_t & mesh, double *buffer, double mask, double t, double p, int m, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_intpl_LGP1_template(mesh,buffer,mask,t,p,m,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_intpl_LGP1(mesh_t & mesh, fcomplex *buffer, fcomplex mask, double t, double p, int m, fcomplex *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_intpl_LGP1_template(mesh,buffer,mask,t,p,m,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T>  int fe_intpl_NCP1_template(const mesh_t & mesh,const T *buffer,T mask,double t,double p,int m,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T zz;
  double x, y, beta[3];
  int nodes[3];
  int i;
  int status=0;

  status=fe_test_area(mesh,m);
  if(status)TRAP_ERR_RETURN(status,1,"%s:fe_test_area(,%d)=%d\n",__func__,m,status);

//  status = fe_affine_directe(t, p, &x, &y, mesh.type);
  status = fe_affine_directe(mesh.triangles[m],t, p, &x, &y, mesh.type);
  if(status != 0)
    goto error;

  fe_NCP1base(x, y, beta);
  for(i = 0; i < 3; i++) {
    nodes[i] = mesh.triangles[m].edges[i];
    }

  zz = 0;
  for(i = 0; i < 3; i++) {
    zz += buffer[nodes[i]] * (T) beta[i];
    }

  *z = zz;
  return (status);

error:
  printf("error in fe_intpl_NCP1() \n");
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_intpl_DGP1_template(const mesh_t & mesh,const T *buffer,T mask,double t,double p,int m,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T zz;
  double x, y, beta[3];
  int i;
  int status;

  status=fe_test_area(mesh,m);
  if(status)TRAP_ERR_RETURN(status,1,"%s:fe_test_area(,%d)=%d\n",__func__,m,status);

  status = fe_affine_directe(mesh.triangles[m],t, p, &x, &y, mesh.type);
  if(status != 0)
    goto error;

  status = fe_LGPbase(mesh, LGP1, x, y, beta);

  zz = 0;
  for(i = 0; i < 3; i++) {
    zz += buffer[m*3+i] * (T) beta[i];
    }

  *z = zz;
  return (status);

error:
  printf("error in fe_intpl_D_LGP1() \n");
  return (-1);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   template <typename T> int fe_intpl_NCP1_template(mesh_t mesh, T *buffer, double t, double p, int m, T *z)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   T zz;
//   double x, y, beta[3];
//   int nodes[3];
//   int i, j, k, node;
//   int status=0;
//   int nndes, elt;
//
//   status = fe_initaffine(&mesh, m);
//   if(status != 0)
//     goto error;
//
// //  status = fe_affine_directe(t, p, &x, &y, mesh.type);
//   status = fe_affine_directe(mesh.triangles[m],t, p, &x, &y, mesh.type);
//   if(status != 0)
//     goto error;
//
//   fe_NCP1base(x, y, beta);
//   for(i = 0; i < 3; i++) {
//     nodes[i] = mesh.triangles[m].edges[i];
//     }
//
//   zz = 0;
//   for(i = 0; i < 3; i++) {
//     zz += buffer[nodes[i]] * beta[i];
//     }
//
//   *z = zz;
//   return (status);
//
// error:
//   printf("error in fe_intpl_NCP1() \n");
//   return (-1);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_intpl_DNP1_template(const mesh_t & mesh,const T *buffer,T mask,double t,double p,int m,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T zz;
  double x, y, beta[3];
  int i;
  int status;

  status=fe_test_area(mesh,m);
  if(status)TRAP_ERR_RETURN(status,1,"%s:fe_test_area(,%d)=%d\n",__func__,m,status);

  status = fe_affine_directe(mesh.triangles[m],t, p, &x, &y, mesh.type);
  if(status != 0)
    goto error;

  fe_NCP1base(x, y, beta);

  zz = 0;
  for(i = 0; i < 3; i++) {
    zz += buffer[m*3+i] * (T) beta[i];
    }

  *z = zz;
  return (status);

error:
  printf("error in fe_intpl_D_NCP1() \n");
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_intpl_LGP2_template(const mesh_t & mesh,const T *buffer,T mask,double t,double p,const int m,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T zz;
  double x, y, beta[6];
  int nodes[6];
  int i;
  int status;

  status=fe_test_area(mesh,m);
  if(status) TRAP_ERR_RETURN(status,1,"%s:fe_test_area(,%d)=%d\n",__func__,m,status);

  status = fe_affine_directe(mesh.triangles[m],t, p, &x, &y, mesh.type);
  if(status) TRAP_ERR_RETURN(status,1,"%s:fe_affine_directe() error\n",__func__,m);

  fe_LGP2base(x, y, beta);
  for(i = 0; i < 6; i++) {
    nodes[i] = mesh.LGP2descriptor.NIbE[m][i];
    }

  zz = 0;
  for(i = 0; i < 6; i++) {
    zz += buffer[nodes[i]] * (T) beta[i];
    }

//  zz = buffer[nodes[0]];
  *z = zz;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_intpl2D(const mesh_t & mesh,discretisation_t descriptor,T *buffer,double t,double p,int m,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T zz;
  double x, y, beta[descriptor.nnpe];
  int nodes[6];
  int i, j, k, node;
  int status;
  int nndes, elt;

  status=fe_test_area(mesh,m);
  if(status)TRAP_ERR_RETURN(status,1,"%s:fe_test_area(,%d)=%d\n",__func__,m,status);

  status = fe_affine_directe(mesh.triangles[m],t, p, &x, &y, mesh.type);
  if(status != 0)
    goto error;

  fe_LGPbase(mesh, descriptor.type, x, y, beta);
  for(i = 0; i < 6; i++) {
    nodes[i] = descriptor.NIbE[m][i];
    }

  zz = 0;
  for(i = 0; i < 6; i++) {
    zz += buffer[nodes[i]] * beta[i];
    }

  *z = zz;
  return (status);

error:
  printf("error in fe_intpl2D() \n");
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_interpolate2D_template(const mesh_t & mesh,int discretisation,const T *buffer,T mask,double t,double p,const int m,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  if(m<0){
    *z=mask;
    return EINVAL;
    }
  
  switch(discretisation) {
    case LGP0:
      status = fe_intpl_LGP0_template(mesh, buffer, mask, t, p, m, z);
      break;

    case LGP1:
      status = fe_intpl_LGP1_template(mesh, buffer, mask, t, p, m, z);
      break;

    case DGP1:
      status = fe_intpl_DGP1_template(mesh, buffer, mask, t, p, m, z);
      break;

    case LGP2:
      status = fe_intpl_LGP2_template(mesh, buffer, mask, t, p, m, z);
      break;

    case NCP1:
      status = fe_intpl_NCP1_template(mesh, buffer, mask, t, p, m, z);
      break;

    case DNP1:
      status = fe_intpl_DNP1_template(mesh, buffer, mask, t, p, m, z);
      break;

    case QLP1:
      status = fe_intpl_QLP1(mesh, buffer, mask, t, p, m, z);
      break;

    case CQP0:
//      status = fe_intpl_CQP0(mesh, buffer, mask, t, p, m, z);
      *z=buffer[m];
      status=0;
      break;

    case CQN1:
      status = fe_intpl_CQN1(mesh, buffer, mask, t, p, m, z);
      break;

    default:
      TRAP_ERR_EXIT(ENOEXEC,"%d discretisation not known\n",discretisation);

    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_interpolate2D(const mesh_t & mesh,int discretisation,const float *buffer,float mask,double t,double p,const int m,float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_interpolate2D_template(mesh,discretisation,buffer,mask,t,p,m,z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_interpolate2D(const mesh_t & mesh,int discretisation,const double *buffer,double mask,double t,double p,const int m,double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_interpolate2D_template(mesh,discretisation,buffer,mask,t,p,m,z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_interpolate2D(const mesh_t & mesh,int discretisation,const complex<float> *buffer,complex<float> mask,double t,double p,const int m,complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_interpolate2D_template(mesh,discretisation,buffer,mask,t,p,m,z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_interpolate2D(const mesh_t & mesh,int discretisation,const double *buffer,double mask,double t,double p,const int m,complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double d;
  status=fe_interpolate2D_template(mesh,discretisation,buffer,mask,t,p,m,&d);
  *z=d;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_interpolate2D(const mesh_t & mesh,int discretisation,const complex<double> *buffer,complex<double> mask,double t,double p,const int m,complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_interpolate2D_template(mesh,discretisation,buffer,mask,t,p,m,z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> double fe_extrapolate_template(const discretisation_t &descriptor,const complex<T> *tide,double x,double y,double dmax,complex<T> *zz,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// try to use nearest node in mesh
/**
\param dmax maximum distance in m
\param *zz pointer to complex amplitude that will be modified only if the distance is whithin \c dmax
\return distance to nearest node
*/
/*----------------------------------------------------------------------------*/
{
  const int node=fe_nearest_node (descriptor, x, y, 0, 0);
  const double d=fe_distance(descriptor,node,x, y);
  
  if(d<dmax){
    if(verbose==1) printf("use nearest node %d: lon=%9.3lf lat=%9.3lf\n", node, descriptor.nodes[node].lon, descriptor.nodes[node].lat);
    *zz=tide[node];
    }
  
  return d;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_extrapolate(const discretisation_t &descriptor,const complex<float> *tide,double x,double y,double dmax,complex<float> *zz,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_extrapolate_template(descriptor,tide,x,y,dmax,zz,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fe_extrapolate(const discretisation_t &descriptor,const complex<double> *tide,double x,double y,double dmax,complex<double> *zz,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_extrapolate_template(descriptor,tide,x,y,dmax,zz,verbose);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> int fe_mapVdepth_template(mesh_t & mesh, grid_t grid,  T *c, T mask, int algorithm, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  interpolate unstructured depths on grid
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  discretisation_t descriptor;
  int   *elements;
  float *buffer;

  for(size_t j = 0; j < grid.ny; j++) {
    for(size_t i = 0; i < grid.nx; i++) {
      size_t m = i + j * grid.nx;
      c[m] = mask;
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  identify elements corresponding to grid nodes

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(verbose==1) printf("identify elements, %d elements, %d structured grid nodes\n",mesh.ntriangles,grid.Hsize());
  
  elements=fe_scan_elements(mesh,grid,algorithm);
  
  if(verbose==1) printf("identify elements done\n");
  
  size_t count=occurence(-1, elements, grid.Hsize());
  printf("%s : %d mesh-inside points over %d\n", __func__, grid.Hsize()-count,grid.Hsize());
      
  buffer=new float[mesh.nvtxs];
  
  for(size_t n=0;n<mesh.nvtxs;n++) buffer[n]=mesh.vertices[n].h;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  interpolate

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(verbose==1) printf("interpolate depths\n");
#pragma omp parallel for
  for(size_t j = 0; j < grid.ny; j++) {
    for(size_t i = 0; i < grid.nx; i++) {
      size_t m = i + j * grid.nx;
      if(elements[m] != -1) {
        double t, p;
        float z;
        grid.xy(m, t, p);
        status=fe_intpl_LGP1(mesh, buffer, mask, t, p, elements[m], &z);
        c[m]=-z;
        }
      else {
        c[m] = mask;
        }
      }
    }
  
  delete[] buffer;
  delete[] elements;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_mapVdepth(mesh_t & mesh, grid_t grid,  float *c, float mask, int algorithm, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=fe_mapVdepth_template( mesh, grid,  c, mask, algorithm, verbose);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_mapVdepth(mesh_t & mesh, grid_t grid,  short *c, short mask, int algorithm, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=fe_mapVdepth_template( mesh, grid,  c, mask, algorithm, verbose);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename B,typename T> int fe_map_template(const mesh_t & mesh,const B *buffer,int discretisation,const grid_t & grid,int *elts,T *c,T mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  interpolate unstructured buffer on structured grid
 
  mesh             : unstructured grid
  *buffer          : unstructured data
  discretisation   : discretisation
  grid             : structured grid
  *elts            : output of fe_scan_elements(), generally fe_scan_elements(mesh,grid,0)
                     array[grid.ny*grid.nx] of element indexes, or -1 if outside
  *c               : output buffer for the structured data. It must be allocated
                     to the number of nodes in the structured grid!
  mask             : mask value

  returns always 0
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int j;

  for(j = 0; j < grid.ny; j++) {
    for(int i = 0; i < grid.nx; i++) {
      int m = i + j * grid.nx;
      c[m] = mask;
      }
    }
  
  int nprocs=initialize_OPENMP(-1, 0);
#pragma omp parallel for if(nprocs>1)
  for(j = 0; j < grid.ny; j++) {
    for(int i = 0; i < grid.nx; i++) {
      double t, p;
      int m = i + j * grid.nx;
//       t = grid.x[m];
//       p = grid.y[m];
      grid.xy(i,j,t,p);
      status=fe_interpolate2D(mesh, discretisation, buffer, mask, t, p, elts[m], &(c[m]));
      }
    }

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_map(const mesh_t & mesh,float *buffer,int discretisation,const grid_t & grid,int *elts,float *c,float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_map_template(mesh,buffer,discretisation,grid,elts,c,mask);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_map(const mesh_t & mesh,double *buffer,int discretisation,const grid_t & grid,int *elts,double *c,double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_map_template(mesh,buffer,discretisation,grid,elts,c,mask);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_map(const mesh_t & mesh,complex<float> *buffer,int discretisation,const grid_t & grid,int *elts,complex<float> *c,complex<float> mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_map_template(mesh,buffer,discretisation,grid,elts,c,mask);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_map(const mesh_t & mesh, complex<double> *buffer, int discretisation, grid_t grid, int *elts, complex<double> *c, complex<double> mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_map_template(mesh,buffer,discretisation,grid,elts,c,mask);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> int tide_atlasUG2positions_template(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, complex<T> *z,complex<T> cmask,const mesh_t &mesh,int discretisation,int verbose=1,const char *wave=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///extract complex amplitudes of unstructured atlas at given locations
/**
\param verbose Default:1
\sa tide_atlas2positions_template() tide_atlasSG2positions_template()
*/
/*----------------------------------------------------------------------------*/
{
  int i, status;
  const discretisation_t *descriptor;
  int frame=-1;
  complex<T> *tide=NULL;
  int *elements;
  /*
    + takes ages for predictor
    + does not take ages for 
  */
#define fe_extrapolate_TAKES_AGES 1
#if !fe_extrapolate_TAKES_AGES
  const double dmax=10e3;/*< maximum distance in m */
#endif
  int count,count2;
  
  descriptor=get_descriptor_address(mesh,discretisation);
  
  exitIfNull(tide=new complex<T>[descriptor->nnodes]);
  
  if(wave!=NULL){
    status=poc_get_frame_from_name(filename,wave,&frame,verbose);
    }
  
  status=poc_get_cvara(filename,v1,v2,frame,tide,verbose);
  if(status!=NC_NOERR){
    aset(z,npositions,cmask);
    NC_TRAP_ERROR(return,status,verbose,"poc_get_cvara(\"%s\",\"%s\",\"%s\",%d,) error",filename,v1,v2,frame);
    }
  
  elements=fe_detect(mesh,lon,lat,npositions);
  
  int nnegative=0;
  for(i=0;i<npositions;i++)
    if(elements[i]<0)
      nnegative++;
  
  count=0;
  count2=0;
  double d=+INFINITY,lon0=lon[0],lat0=lat[0];
  
  //#pragma omp parallel for private(i,d,lon0,lat0,status,count2) reduction(+:count)
  for(i=0;i<npositions;i++){
    const int *eli=&elements[i];
    complex<T> *zi=&z[i];
    const double *loni=&lon[i],*lati=&lat[i];
    
    if(*eli<0){
      count2++;
#if !fe_extrapolate_TAKES_AGES
fprintf(stderr,el);STDERR_BASE_LINE("%9g%% (%9g;%9g)\r",count2*100./nnegative,*loni,*lati);
      /*
      if(d>dmax)
        d=geo_distance(lon0,lat0,*loni,*lati,'m');
      if(d<dmax)
        d=fe_extrapolate(descriptor,tide,*loni,*lati,dmax,zi,0);
      */
      d=fe_extrapolate(*descriptor,tide,*loni,*lati,dmax,zi,0);
      if(d>dmax){
        *zi=cmask;
        //TRAP_ERR_EXIT(ENOEXEC,"testing\n");
        }
#else
      *zi=cmask;
#endif
      status=0;
      }
    else{
      status=fe_interpolate2D(mesh,discretisation,tide,cmask,*loni,*lati,*eli,zi);
      d=0;
      }
    
    lon0=*loni;
    lat0=*lati;
    
    if(verbose>0 && status)count++;
    }
  
  if(verbose>0 && count>0)STDERR_BASE_LINE("fe_interpolate2D() error on %d/%d points\n",count,npositions);
  
  status=NC_NOERR;
  
  delete[]elements;
  delete[]tide;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_atlasUG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, complex<double> *z,complex<double> cmask,const mesh_t &mesh,int discretisation,int verbose,const char *wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlasUG2positions_template(filename,v1,v2,lon,lat,npositions,z,cmask,mesh,discretisation,verbose,wave);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_atlasUG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, complex<float> *z,complex<float> cmask,const mesh_t &mesh,int discretisation,int verbose,const char *wave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlasUG2positions_template(filename,v1,v2,lon,lat,npositions,z,cmask,mesh,discretisation,verbose,wave);
  return status;
}
