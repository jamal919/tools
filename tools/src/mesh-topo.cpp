
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

*******************************************************************************/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "archive.h"
#include "map.h"
#include "topo.h"
#include "polygones.h"
#include "grd.h"
#include "functions.h"
#include "discretisation.h"
#include "mass.h"
#include "matrix.h"
#include "netcdf-classes.h"
#include "netcdf-proto.h"
#include "datastream.h"

#include "map.def"

#include "zapper.h"

#include "version-macros.def" //for VERSION and REVISION

extern void fe_LGP0_to_LGP1(double *in, double *out, mesh_t mesh);

extern  discretisation_t initCQP1(mesh_t mesh, int type);
extern  discretisation_t initCQP0(mesh_t mesh, int type);

int field_interpolation(SGfield_t<float> & field, double x, double y, float *z);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

 \param prog_name name of the program to be printed by this help function : as[0]
 */
/*----------------------------------------------------------------------------*/
{
    /** The code of the body of this function is :
     \code /**/ // COMPILED CODE BELOW !!!
    printf("\n"
            " NAME AND VERSION\n"
            "  %s version " VERSION " " REVISION "\n", prog_name);
    printf("\n USE");
    printf("\n   %s -m MESH_FILE -b BATHYMETRY_FILE -o OUTPUT_FILE -d DISCRETISATION -p PAIRE  \n",prog_name);
            printf("\n DESCRIPTION");
            printf(
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
                    "\n   Interpolate a bathymetry and its gradient on a given FE mesh");
            printf("\n ARGUMENTS :\n"
                    "   -m       the mesh on which to interpolate the bathymetry \n"
                    "            (.nei format, if needed, use mesh-format to convert from gmesh or other)\n"
/** ***************************************************
 \todo Jan 8, 2013 :  Clement MAYET : precise the netcdf format needed
*************************************************** */
                    "   -b       bathymetry file, must be a netcdf or grd file\n"
                    "   -o       **NOT USED** output file name \n"
/** ***************************************************
 \todo Jan 8, 2013 :  Clement MAYET : clarify the following option
*************************************************** */
                    "   -d       **NOT USED** discretisation \n"
                    "   -p       computational element paire (LGP0xLGP1 DGP1xLGP2 DNP1xLGP2 Q1xQ0)\n"
                    "            used to compute optimal topography, the topo will be saved on elevation nodes \n "
                    "            and the topography gradient on velocity nodes \n"
                    "   -h       print this help\n");

            TRAP_ERR_EXIT(-1,"exiting\n"); /** \endcode */
        }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_edgetable_quick(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,k,m,n1,n2,n;
  int   nedges,count=0;
  int   nelt,nndes;
  int   *connexions;
  double t1,t2,p1,p2;
  edge_t *edges;
  triangle_t *elt;

  elt=mesh->triangles;
  nelt=mesh->ntriangles;
  nndes=mesh->nvtxs;

  nedges=0;

/*-----------------------------------------------------------------------------
  count edges*/
  for(n=0; n<nndes; n++) {
    for(j=0; j<mesh->vertices[n].nngh;j++) {
      n2=mesh->vertices[n].ngh[j];
      if(n2>n) nedges++;
      }
    }

  edges=new edge_t[nedges];

  connexions=new int[mesh->nvtxs*mesh->nvtxs];

  for(n1=0; n1<nndes; n1++) {
    for(j=0; j<mesh->vertices[n1].nngh;j++) {
      n2=mesh->vertices[n1].ngh[j];
      if(n2>n1) {
/*-----------------------------------------------------------------------------
        create a new edge entry*/
        edges[count].extremity[0]=n1;
        edges[count].extremity[1]=n2;
        connexions[n2*mesh->nvtxs+n1]=count;
        connexions[n1*mesh->nvtxs+n2]=count;
        edges[count].nshared=0;
        t1=mesh->vertices[n1].lon;
        t2=mesh->vertices[n2].lon;
        if(mesh->type==0) t2=geo_recale(t2,t1,180.0);
        p1=mesh->vertices[n1].lat;
        p2=mesh->vertices[n2].lat;
        edges[count].lon=0.5*(t1+t2);
        edges[count].lat=0.5*(p1+p2);
/*-----------------------------------------------------------------------------
        scan elements shared by n1 and n2*/
        count++;
        }
      }
    }
  nedges=count;

  for(m=0;m<mesh->ntriangles;m++){
    for(i=0;i<3;i++) {
      n1=mesh->triangles[m].vertex[i];
      j=(i+1)%3;
      n2=mesh->triangles[m].vertex[j];
      n=connexions[n2*mesh->nvtxs+n1];
      edges[n].shared[edges[n].nshared]=m;
      edges[n].nshared++;
      k=(i+2)%3;
      mesh->triangles[m].edges[k]=n;
      }
    }

  mesh->nedges=nedges;
  mesh->edges=edges;

  delete[] connexions;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_double_quick(mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,n1,n2,n,status;
  int   nedges;
  int   nelts,nndes;
  edge_t *edges;
  triangle_t *elt;
  double /*t1,t2,p1,p2,*/h1,h2;
  mesh_t work;

  elt=mesh.triangles;

  nelts =mesh.ntriangles;
  nndes =mesh.nvtxs;
  nedges=mesh.nedges;
  edges =mesh.edges;

  work.type=0;

/*-----------------------------------------------------------------------------
  final mesh : original vertices + 1 per edges vertex */
  work.vertices=new vertex_t[nedges+nndes];

  for (n=0;n<nndes;n++) {
    work.vertices[n].lon=mesh.vertices[n].lon;
    work.vertices[n].lat=mesh.vertices[n].lat;
    work.vertices[n].h=mesh.vertices[n].h;
    work.vertices[n].nngh=0;
    }

  for (n=0;n<nedges;n++) {
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
/*-----------------------------------------------------------------------------
    create a new node at the middle of the edge */
//     t1=mesh.vertices[n1].lon;
//     t2=mesh.vertices[n2].lon;
//     if(mesh.type==0) t2=geo_recale(t2,t1,180.0);
//     p1=mesh.vertices[n1].lat;
//     p2=mesh.vertices[n2].lat;
    h1=mesh.vertices[n1].h;
    h2=mesh.vertices[n2].h;
//     work.vertices[n+nndes].lon=0.5*(t1+t2);
//     work.vertices[n+nndes].lat=0.5*(p1+p2);
    work.vertices[n+nndes].lon=edges[n].lon;
    work.vertices[n+nndes].lat=edges[n].lat;
    work.vertices[n+nndes].h=0.5*(h1+h2);
    work.vertices[n+nndes].nngh=0;
    }

 work.nvtxs=nedges+nndes;
 work.ntriangles=0;
 work.triangles= new triangle_t[4*nelts];

/*-----------------------------------------------------------------------------
  list old element to create new elements */
  for(k=0; k<nelts; k++) {
    work.triangles[work.ntriangles].vertex[0]=elt[k].vertex[0];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[2];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[1];

    work.ntriangles++;
    work.triangles[work.ntriangles].vertex[0]=elt[k].vertex[1];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[0];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[2];

    work.ntriangles++;
    work.triangles[work.ntriangles].vertex[0]=elt[k].vertex[2];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[1];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[0];

    work.ntriangles++;
    work.triangles[work.ntriangles].vertex[0]=nndes+elt[k].edges[2];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[0];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[1];

    work.ntriangles++;
    }

/*-----------------------------------------------------------------------------
  rebuild neighbour list */
  status=fe_e2n(&work);
  status=fe_edgetable_quick(&work);

  return(work);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  mesh_t fe_splitelement(mesh_t mesh, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,k,m,n,status;
  int   count=0,n1,n2;
  double t1,p1,t2,p2;
  mesh_t work;

  work.type=mesh.type;

  work.triangles=new triangle_t[4];
  work.vertices=new vertex_t[6];

  for (k=0;k<3;k++) {
    n=mesh.triangles[target].vertex[k];
    work.vertices[count].lon=mesh.vertices[n].lon;
    work.vertices[count].lat=mesh.vertices[n].lat;
    work.vertices[count].h=mesh.vertices[n].h;
    work.vertices[count].code=0;
    work.vertices[count].nngh=0;
    count++;
    }

  for (k=0;k<3;k++) {
    n=mesh.triangles[target].edges[k];
    n1=mesh.edges[n].extremity[0];
    n2=mesh.edges[n].extremity[1];
    t1=mesh.vertices[n1].lon;
    t2=mesh.vertices[n2].lon;
    p1=mesh.vertices[n1].lat;
    p2=mesh.vertices[n2].lat;
    t2=geo_recale(t2,t1,180.0);
    work.vertices[count].lon=0.5*(t1+t2);
    work.vertices[count].lat=0.5*(p1+p2);
    work.vertices[count].h=0.5*(mesh.vertices[n1].h+mesh.vertices[n2].h);
    work.vertices[count].code=0;
    work.vertices[count].nngh=0;
    count++;
    }

  m=0;
  work.triangles[m].vertex[0]=0;
  work.triangles[m].vertex[1]=4;
  work.triangles[m].vertex[2]=5;

  m++;
  work.triangles[m].vertex[0]=5;
  work.triangles[m].vertex[1]=1;
  work.triangles[m].vertex[2]=3;

  m++;
  work.triangles[m].vertex[0]=3;
  work.triangles[m].vertex[1]=2;
  work.triangles[m].vertex[2]=4;

  m++;
  work.triangles[m].vertex[0]=3;
  work.triangles[m].vertex[1]=4;
  work.triangles[m].vertex[2]=5;

  work.ntriangles=4;
  work.nvtxs=6;

  status= fe_e2n(&work);

  for (m=0;m<work.ntriangles;m++) {
    for (i=0;i<3;i++) {
      n=work.triangles[m].vertex[i];
      }
    fe_initaffine(&work,m);
    }

  return(work);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_average(mesh_t mesh, int target, grid_t topogrid, float **buffers, int nbuffers,float mask,double r, double *sum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,k,m,n,status;
  int   kmax=3;
  double t,p,rr,rate,w;
  float h;
  mesh_t work,refined[10];

  rr=0;
  for (i=0;i<3;i++) {
    n=mesh.triangles[target].edges[i];
    rr=max(rr,mesh.edges[n].L);
    }

  rate=rr/r;
  kmax=0;
  while(rate > 0.5) {
    rate=rate/2.0;
    kmax++;
    }

  refined[0]=fe_splitelement(mesh, target);
  status= fe_edgetable(&(refined[0]),0,0);
  for (k=0;k<kmax;k++) {
    rr=0;
    refined[k+1]=fe_double_quick(refined[k]);
    (refined[k]).destroy();
    }

  work=refined[kmax];
  for (m=0;m<work.ntriangles;m++) {
    for (i=0;i<3;i++) {
      n=work.triangles[m].vertex[i];
      }
    fe_initaffine(&work,m);
    }

  for(k=0;k<nbuffers;k++) sum[k]=0.;
  w=0.;
  for (m=0;m<work.ntriangles;m++) {
    t=0;
    p=0;
    for (i=0;i<3;i++) {
      n=work.triangles[m].vertex[i];
      t+=work.vertices[n].lon/3.0;
      p+=work.vertices[n].lat/3.0;
      }
    w+=work.triangles[m].Area;
    for(k=0;k<nbuffers;k++) {
      status=map_interpolation(topogrid,buffers[k],mask,t,p,&h);
      sum[k]+=h*work.triangles[m].Area;
      if(isnan(sum[k])==1) {
        sum[k]=0.0;
        }
      }
    }

  for(k=0;k<nbuffers;k++) sum[k]/=mesh.triangles[target].Area;

  (refined[kmax]).destroy();

  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_averageQ(mesh_t mesh, int target, grid_t topogrid, float **buffers, int nbuffers,float mask,double r, double *sum)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,k,l,/*m,*/n,status;
  int   kmax=3;
  double t,p,rr,rate;
  float h;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;

  gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);

  rr=0;
  for (i=0;i<3;i++) {
    n=mesh.triangles[target].edges[i];
    rr=max(rr,mesh.edges[n].L);
    }

  rate=rr/r;
  kmax=0;
  while(rate > 0.5) {
    rate=rate/2.0;
    kmax++;
    }

  for(l=0;l<nbuffers;l++) sum[l]=0.;
  for (k=0;k<gauss_n;k++) {
//     t=0;
//     p=0;
//     for (i=0;i<3;i++) {
//       n=work.triangles[m].vertex[i];
//       t+=work.vertices[n].lon/3.0;
//       p+=work.vertices[n].lat/3.0;
//       }
    status=fe_affine_inverse(mesh.triangles[target], &t, &p, gauss_x[k], gauss_y[k]);
    t=map_recale(topogrid,t);
    for(l=0;l<nbuffers;l++) {
      status=map_interpolation(topogrid,buffers[l],mask,t,p,&h);
      if(h==mask) {
        printf("fe_averageQ: trouble\n");
        }
      sum[l]+=h*gauss_w[k]*2.0;
      }
    }

  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RHS_generic_T_obsolete(mesh_t mesh, discretisation_t descriptor, double *rhs, grid_t topogrid, float *buffers,float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Triangle Gauss integration

----------------------------------------------------------------------*/
  int i, k, m, n, status;
  double C, area;
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double **base;
  float h;
  int discretisation=descriptor.type;
  
  status=gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);
  
  for(n = 0; n < descriptor.nnodes; n++)
    rhs[n] = 0;
  
  base=new double*[gauss_n];
  for (k=0;k<gauss_n;k++) {
    base[k]=new double[descriptor.nnpe];
    status=fe_LGPbase(mesh, discretisation, gauss_x[k], gauss_y[k], base[k]);
    }

  for(m = 0; m < mesh.ntriangles; m++) {
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n = descriptor.NIbE[m][i];
      for (k=0;k<gauss_n;k++) {
        status=fe_affine_inverse(mesh.triangles[m], &t, &p, gauss_x[k], gauss_y[k]);
        t=map_recale(topogrid,t);
        status=map_interpolation(topogrid,buffers,mask,t,p,&h);
//        h=-1000;
        if(h==mask) {
          printf("RHS_generic_T_obsolete : trouble\n");
          }
        rhs[n] +=h*area*C*gauss_w[k]*base[k][i]*2.0;
        }
      }
    }

  for (k=0;k<gauss_n;k++) {
    delete[] base[k];
    }
  delete[] base;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RHS_generic_Q_obsolete(mesh_t mesh, discretisation_t descriptor, double *rhs, grid_t topogrid, float *buffers,float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Triangle Gauss integration

----------------------------------------------------------------------*/
  int i, k, m, n, status;
  double C, area;
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double **base;
  float h;
  int discretisation=descriptor.type;
  extern int fe_affine_inverse(mesh_t & mesh, quadrangle_t e, double *t, double *p,double x,double y);
  
  status=gauss_init_Q(gauss_n,gauss_x,gauss_y,gauss_w);
  
  for(n = 0; n < descriptor.nnodes; n++)
    rhs[n] = 0;
  
  base=new double*[gauss_n];
  for (k=0;k<gauss_n;k++) {
    base[k]=new double[descriptor.nnpe];
    status=fe_LGPbase(mesh, discretisation, gauss_x[k], gauss_y[k], base[k]);
    }

  for(m = 0; m < mesh.nquadrangles; m++) {
    C=mesh.quadrangles[m].cosinus;
    area = mesh.quadrangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n = descriptor.NIbE[m][i];
      for (k=0;k<gauss_n;k++) {
        status=fe_affine_inverse(mesh, mesh.quadrangles[m], &t, &p, gauss_x[k], gauss_y[k]);
        t=map_recale(topogrid,t);
        status=map_interpolation(topogrid,buffers,mask,t,p,&h);
//        h=-1000;
        if(h==mask) {
          printf("RHS_generic_Q_obsolete : trouble\n");
          }
        rhs[n] +=h*area*C*gauss_w[k]*base[k][i];
        }
      }
    }

  for (k=0;k<gauss_n;k++) {
    delete[] base[k];
    }
  delete[] base;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RHS_generic_obsolete(mesh_t mesh, discretisation_t descriptor, double *rhs, grid_t topogrid, float *buffers,float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
  switch(mesh.ntriangles) {
    case 0:
      RHS_generic_Q_obsolete(mesh, descriptor, rhs, topogrid, buffers, mask);
      break;
    default:
      RHS_generic_T_obsolete(mesh, descriptor, rhs, topogrid, buffers, mask);
      break;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RHS_generic_T(mesh_t & mesh, discretisation_t & descriptor, double *rhs, Scombo_t combo, float missing)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Triangle Gauss integration

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, k, l, m, n, status;
  double C, area;
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double **base;
  float h,hh;
  int discretisation=descriptor.type;
  float mask;
  
  status=gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);
  
  for(n = 0; n < descriptor.nnodes; n++)
    rhs[n] = 0;
  
  base=new double*[gauss_n];
  for (k=0;k<gauss_n;k++) {
    base[k]=new double[descriptor.nnpe];
    status=fe_LGPbase(mesh, discretisation, gauss_x[k], gauss_y[k], base[k]);
    }

//   mask=1.e+10;
  mask=missing;
  
  for(m = 0; m < mesh.ntriangles; m++) {
    vector<int> masked;
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n = descriptor.NIbE[m][i];
      for (k=0;k<gauss_n;k++) {
        status=fe_affine_inverse(mesh.triangles[m], &t, &p, gauss_x[k], gauss_y[k]);
        h=mask;
        for(l=0;l<combo.nfields;l++) {
          t=map_recale(*combo.fields[l].grid,t);
          status=field_interpolation(combo.fields[l],t,p,&hh);
          if(hh==combo.fields[l].mask) {
            status=map_interpolation(*combo.fields[l].grid, combo.fields[l].x, combo.fields[l].mask,t,p,&hh);
            }
          if(hh!=combo.fields[l].mask) {
/**----------------------------------------------------------------------------
            give priority to 1st fields */
            if(h==mask) h=hh;
            }
          }
        if(h==mask) {
          masked.push_back(k);
          h=0;
          }
        rhs[n] +=h*area*C*gauss_w[k]*base[k][i]*2.0;
        }
      }
    if(masked.size()!=0) {
      printf("RHS_generic_T: masked value at element %d lon=%lf lat=%lf (set to %f)\n",m,t,p,missing);
      }
    }

  for (k=0;k<gauss_n;k++) {
    delete[] base[k];
    }
  delete[] base;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RHS_generic_Q(mesh_t & mesh, discretisation_t & descriptor, double *rhs, Scombo_t combo, float missing)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Quadrangle Gauss integration

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, k, l=0, m, n, status;
  double C, area;
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double **base;
  float h;
  int discretisation=descriptor.type;
  extern int fe_affine_inverse(mesh_t & mesh, quadrangle_t e, double *t, double *p,double x,double y);
  
  status=gauss_init_Q(gauss_n,gauss_x,gauss_y,gauss_w);
  
  for(n = 0; n < descriptor.nnodes; n++)
    rhs[n] = 0;
  
  base=new double*[gauss_n];
  for (k=0;k<gauss_n;k++) {
    base[k]=new double[descriptor.nnpe];
    status=fe_LGPbase(mesh, discretisation, gauss_x[k], gauss_y[k], base[k]);
    }

  for(m = 0; m < mesh.nquadrangles; m++) {
    C=mesh.quadrangles[m].cosinus;
    area = mesh.quadrangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n = descriptor.NIbE[m][i];
      for (k=0;k<gauss_n;k++) {
        status=fe_affine_inverse(mesh, mesh.quadrangles[m], &t, &p, gauss_x[k], gauss_y[k]);
        for(k=0;k<combo.nfields;k++) {
          t=map_recale(*combo.fields[k].grid,t);
          status=field_interpolation(combo.fields[k],t,p,&h);
          if(h==combo.fields[l].mask) {
            printf("RHS_generic_Q : trouble\n");
            }
          }
        rhs[n] +=h*area*C*gauss_w[k]*base[k][i];
        }
      }
    }

  for (k=0;k<gauss_n;k++) {
    delete[] base[k];
    }
  delete[] base;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void RHS_generic(mesh_t & mesh, discretisation_t & descriptor, double *rhs, Scombo_t combo, float missing)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
  switch(mesh.ntriangles) {
    case 0:
      RHS_generic_Q(mesh, descriptor, rhs, combo, missing);
      break;
    default:
      RHS_generic_T(mesh, descriptor, rhs, combo, missing);
      break;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *CHK_generic(mesh_t & mesh, discretisation_t & descriptor, double *UGtopo, grid_t topogrid, float *buffers,float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Mass matrix for P2 integrale scalar product

----------------------------------------------------------------------*/
  int i, k, m, n, status;
  double C, area;
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double **base;
  gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);
  float h,*UGh;
  double *rms;
  double global;
  int discretisation=descriptor.type;
    
  UGh=new float[gauss_n];
  base=new double*[gauss_n];
  for (k=0;k<gauss_n;k++) {
    base[k]=new double[descriptor.nnpe];
    status=fe_LGPbase(mesh, discretisation, gauss_x[k], gauss_y[k], base[k]);
    }

  rms=new double[mesh.ntriangles];
  global=0.0;
  
  for(m = 0; m < mesh.ntriangles; m++) {
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    rms[m]=0.0;
    for (k=0;k<gauss_n;k++) {
      UGh[k]=0.0;
      for(i = 0; i < descriptor.nnpe; i++) {
        n = descriptor.NIbE[m][i];
        UGh[k]+=UGtopo[n]*base[k][i];
        }
      status=fe_affine_inverse(mesh.triangles[m], &t, &p, gauss_x[k], gauss_y[k]);
      t=map_recale(topogrid,t);
      status=map_interpolation(topogrid,buffers,mask,t,p,&h);
      if(h==mask) {
        printf("CHK_generic : trouble\n");
        }
      rms[m]+=(UGh[k]-h)*(UGh[k]-h);
      }
    rms[m]=sqrt(rms[m]/gauss_n);
    }
  
  for (k=0;k<gauss_n;k++) {
    delete[] base[k];
    }
  delete[] base;
  
  return(rms);

}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_topo01_generic_obsolete(const mesh_t & mesh, const discretisation_t & descriptor, grid_t topogrid, float *topo, float topomask, char *dataset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,status;
  double t,p;
  float *depths,*slope_x,*slope_y,mask=99999.9;
  float *topo_x,*topo_y;
  size_t size;
  char *comments[2],*filename;
  int discretisation;

  comments[0]=new char[1024];
  comments[1]=new char[1024];
  filename=new char[1024];

  printf("#################################################################\n");
  printf("interpolate depths\n");
  depths=new float[descriptor.nnodes];

  if(topogrid.modeH!=0)
    topogrid.dx=topogrid.x[1]-topogrid.x[0];

  for (n=0;n<descriptor.nnodes;n++) {
    t=descriptor.nodes[n].lon;
    p=descriptor.nodes[n].lat;
    t=map_recale(topogrid,t);
    status=map_interpolation(topogrid,topo,topomask,t,p,&(depths[n]));
    if(depths[n]==topomask) {
      printf("fe_topo01_generic_obsolete : trouble\n");
      status=map_interpolation(topogrid,topo,topomask,t,p,&(depths[n]));
      }
/**----------------------------------------------------------------------------
    model needs positive depths */
    if(depths[n]!=topomask) {
      depths[n]=-depths[n];
      }
    else {
      depths[n]=mask;
      }
    }
  
  discretisation=descriptor.type;
  
  sprintf(comments[0],"%s bathymetry, unit in meters",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, direct interpolation from %s",dataset);
  sprintf(filename,"topo-%s-0.s2r",discretisation_name(discretisation));
  status=quoddy_saver1(filename, descriptor.nnodes, depths,comments);

//   rms=CHK_generic(mesh, descriptor, rhs, topogrid, topo, topomask);
//   sprintf(varname,"rmsH-%s",discretisation_name(discretisation));
//   status=archiving_UGdummy2D("mesh-topo.check.nc", mesh, varname, rms,(int) LGP0);


  printf("#################################################################\n");
  printf("interpolate depth's gradients\n");
  size=topogrid.nx*topogrid.ny;
  topo_x=new float[size];
  topo_y=new float[size];
  status= map_gradient(topogrid, topogrid.nx, topo, topomask, GEOCENTRIC, topo_x, topo_y);

  slope_x=new float[descriptor.nnodes];
  slope_y=new float[descriptor.nnodes];
  for (n=0;n<descriptor.nnodes;n++) {
    t=descriptor.nodes[n].lon;
    p=descriptor.nodes[n].lat;
    t=map_recale(topogrid,t);
    status=map_interpolation(topogrid,topo_x,topomask,t,p,&(slope_x[n]));
    status=map_interpolation(topogrid,topo_y,topomask,t,p,&(slope_y[n]));
    slope_x[n]*=1000.0;
    slope_y[n]*=1000.0;
    }

  sprintf(comments[0],"%s bathymetry slope, dimensionless",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, direct interpolation from %s",dataset);
  sprintf(filename,"slope-%s-0.v2r",discretisation_name(discretisation));
  status=quoddy_saver2(filename, descriptor.nnodes, slope_x, slope_y,comments);

  delete[] depths;

  delete[] topo_x;
  delete[] topo_y;

  delete[] slope_x;
  delete[] slope_y;

  delete[] comments[0];
  delete[] comments[1];

  printf("%s completed (raw interpolation)\n\n",__func__);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_topo02_generic_obsolete(mesh_t mesh, discretisation_t descriptor, grid_t topogrid, float *topo,float topomask, char *dataset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   n,status;
  float *depths,*slope_x,*slope_y;
  double *rhs;
  float *topo_x,*topo_y;
  size_t size;
  char *comments[2], *filename, *varname;
  mesh_t work;
//   double *dx,*dy,*dz;
  int discretisation;
  double *rms, mask=1.e+10;

  comments[0]=new char[1024];
  comments[1]=new char[1024];
  filename =new char[1024];
  varname  =new char[1024];

  depths=new float[descriptor.nnodes];

//  status=map_resolution( topogrid,  &dx,  &dy,  &dz);
  size=topogrid.nx*topogrid.ny;
  topo_x=new float[size];
  topo_y=new float[size];
  status= map_gradient(topogrid, topogrid.nx, topo, topomask, GEOCENTRIC, topo_x, topo_y);

  slope_x=new float[descriptor.nnodes];
  slope_y=new float[descriptor.nnodes];

  rhs=new double[descriptor.nnodes];
  
  printf("#################################################################\n");
  printf("interpolate depths\n");
  RHS_generic_obsolete(mesh, descriptor, rhs, topogrid, topo, topomask);
  status = LinearSystem_solve(descriptor.massmatrix,rhs);
  for (n=0;n<descriptor.nnodes;n++) {
    depths[n]=-rhs[n];
    }

  discretisation=descriptor.type;

  sprintf(comments[0],"%s bathymetry, unit in meters",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, optimal estimate from %s",dataset);
  sprintf(filename,"topo-%s-1-old.s2r",discretisation_name(discretisation));
  status=quoddy_saver1(filename, descriptor.nnodes, depths,comments);
  
  rms=CHK_generic(mesh, descriptor, rhs, topogrid, topo, topomask);
  sprintf(varname,"rmsH-%s",discretisation_name(discretisation));
  if(mesh.ntriangles!=0) {
    status=archiving_UGdummy2D("mesh-topo.check.nc", mesh, varname, "m", rms, mask, (int) LGP0);
    }
  else {
    status=archiving_UGdummy2D("mesh-topo.check.nc", mesh, varname, "m", rms, mask, (int) CQP0);
    }
  
  printf("#################################################################\n");
  printf("interpolate depth's gradients\n");
  RHS_generic_obsolete(mesh, descriptor, rhs, topogrid, topo_x, topomask);
  status = LinearSystem_solve(descriptor.massmatrix,rhs);
  for (n=0;n<descriptor.nnodes;n++) {
    slope_x[n]=(float) rhs[n];
    slope_x[n]*=1000.0;
    }
  
  RHS_generic_obsolete(mesh, descriptor, rhs, topogrid, topo_y, topomask);
  status = LinearSystem_solve(descriptor.massmatrix,rhs);
  for (n=0;n<descriptor.nnodes;n++) {
    slope_y[n]=(float) rhs[n];
    slope_y[n]*=1000.0;
    }

  sprintf(comments[0],"%s bathymetry slope, dimensionless",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, optimal estimate from %s",dataset);
  sprintf(filename,"slope-%s-1-old.v2r",discretisation_name(discretisation));
  status=quoddy_saver2(filename, descriptor.nnodes, slope_x, slope_y,comments);

  delete[] rhs;

  delete[] depths;

  delete[] topo_x;
  delete[] topo_y;

  delete[] slope_x;
  delete[] slope_y;

  delete[] comments[0];
  delete[] comments[1];

  printf("%s completed (optimal)\n\n",__func__);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_persistence(const mesh_t & mesh, const discretisation_t & descriptor, float *depths, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  size_t missing=occurence(mask, depths, descriptor.nnodes);
  
  while(missing!=0) {
    for (int n=0;n<descriptor.nnodes;n++) {
      if(depths[n]==mask) {
        double sum=0.0, weight=0.0;
        for(int k=0;k<descriptor.nodes[n].nnghbs;k++) {
          int m=descriptor.nodes[n].nghbs[k];
          if(depths[m]!=mask) {
            sum+=depths[m];
            weight+=1.0;
            }
          }
        if(weight!=0) depths[n]=sum/weight;
        }
      }
    missing=occurence(mask, depths, descriptor.nnodes);
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_topo01_generic(const mesh_t & mesh, const discretisation_t & descriptor, Scombo_t combo, char *dataset, float factor, float missing)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,n,status;
  double t,p;
  float *depths,*slope_x,*slope_y,z;
//   float mask=99999.9;
  float mask;
  float *topo_x,*topo_y;
  size_t size;
  char *comments[2],*filename;
  int discretisation;
  bool holes=false;

  comments[0]=new char[1024];
  comments[1]=new char[1024];
  filename=new char[1024];

  printf("#################################################################\n");
  printf("interpolate depths\n");
  mask=missing;
  depths=new float[descriptor.nnodes];
  for (n=0;n<descriptor.nnodes;n++) depths[n]=mask;
  
  for(k=0;k<combo.nfields;k++) {
    combo.fields[k].acquire(0.);
    set_grid_list(combo.fields[k].grid);
    for (n=0;n<descriptor.nnodes;n++) {
      t=descriptor.nodes[n].lon;
      p=descriptor.nodes[n].lat;
      t=map_recale(*combo.fields[k].grid,t);
      status=field_interpolation(combo.fields[k],t,p,&z);
      if(z==combo.fields[k].mask) {
        holes=true;
        }
/**----------------------------------------------------------------------------
      model needs positive depths */
      if(z!=combo.fields[k].mask) {
/**----------------------------------------------------------------------------
        give priority to 1st fields */
        if(depths[n]==mask) depths[n]=factor*z;
        }
      }
    }
  
  status=fe_persistence(mesh, descriptor, depths, mask);
  
  for (n=0;n<descriptor.nnodes;n++) {
    if(depths[n]==mask) {
      t=descriptor.nodes[n].lon;
      p=descriptor.nodes[n].lat;
      printf("%s : hole at node %6d: lon=%lf lat=%lf (set to %f)\n",__func__,n,t,p,missing);
      }
    }
  
  discretisation=descriptor.type;
  
  sprintf(comments[0],"%s bathymetry, unit in meters",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, direct interpolation from %s",dataset);
  sprintf(filename,"topo-%s-0.s2r",discretisation_name(discretisation));
  status=quoddy_saver1(filename, descriptor.nnodes, depths,comments);

//   rms=CHK_generic(mesh, descriptor, rhs, topogrid, topo, topomask);
//   sprintf(varname,"rmsH-%s",discretisation_name(discretisation));
//   status=archiving_UGdummy2D("mesh-topo.check.nc", mesh, varname, rms,(int) LGP0);


  printf("#################################################################\n");
  printf("interpolate depth's gradients\n");
  
//   mask=0.0;
  
  slope_x=new float[descriptor.nnodes];
  slope_y=new float[descriptor.nnodes];
  
  for (n=0;n<descriptor.nnodes;n++) slope_x[n]=mask;
  for (n=0;n<descriptor.nnodes;n++) slope_y[n]=mask;

  for(k=0;k<combo.nfields;k++) {
    size=combo.fields[k].nvalues();
    topo_x=new float[size];
    topo_y=new float[size];
    status= map_gradient(*combo.fields[k].grid,combo.fields[k].grid->nx,combo.fields[k].x,combo.fields[k].mask, GEOCENTRIC, topo_x, topo_y);

    for (n=0;n<descriptor.nnodes;n++) {
      t=descriptor.nodes[n].lon;
      p=descriptor.nodes[n].lat;
      t=map_recale(*combo.fields[k].grid,t);
      status=map_interpolation(*combo.fields[k].grid,topo_x,combo.fields[k].mask,t,p,&z);
      if(z==combo.fields[k].mask) {
        holes=true;
        }
      if(z!=combo.fields[k].mask) {
/**----------------------------------------------------------------------------
        give priority to 1st fields */
        if(slope_x[n]==mask) slope_x[n]=z*1000.0;
        }
      status=map_interpolation(*combo.fields[k].grid,topo_y,combo.fields[k].mask,t,p,&z);
      if(z==combo.fields[k].mask) {
        holes=true;
        }
      if(z!=combo.fields[k].mask) {
/**----------------------------------------------------------------------------
        give priority to 1st fields */
        if(slope_y[n]==mask) slope_y[n]=z*1000.0;
        }
      }
    
    delete[] topo_x;
    delete[] topo_y;
    }
  
  status=fe_persistence(mesh, descriptor, slope_x, mask);
  status=fe_persistence(mesh, descriptor, slope_y, mask);

  for (n=0;n<descriptor.nnodes;n++) {
    if(slope_x[n]==mask or slope_y[n]==mask) {
      t=descriptor.nodes[n].lon;
      p=descriptor.nodes[n].lat;
      printf("%s : hole at node %6d: lon=%lf lat=%lf (set to %f)\n",__func__,n,t,p,mask);
      }
    }
  
  sprintf(comments[0],"%s bathymetry slope, dimensionless",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, direct interpolation from %s",dataset);
  sprintf(filename,"slope-%s-0.v2r",discretisation_name(discretisation));
  status=quoddy_saver2(filename, descriptor.nnodes, slope_x, slope_y,comments);

  for(k=0;k<combo.nfields;k++) {
    combo.fields[k].destroy();
    }

  delete[] depths;

  delete[] slope_x;
  delete[] slope_y;

  delete[] comments[0];
  delete[] comments[1];

  printf("%s completed (raw interpolation)\n",__func__);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_topo02_generic(mesh_t & mesh, discretisation_t & descriptor, Scombo_t & combo, char *dataset, float factor, float missing, char *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,n,status;
  float *depths,*slope_x,*slope_y;
  double *rhs;
  size_t size;
  char *comments[2], *filename, *varname;
  mesh_t work;
//   double *dx,*dy,*dz;
  int discretisation;
//   double *rms;
  Scombo_t combo_gx, combo_gy;

  comments[0]=new char[1024];
  comments[1]=new char[1024];
  filename =new char[1024];
  varname  =new char[1024];

  depths=new float[descriptor.nnodes];

  slope_x=new float[descriptor.nnodes];
  slope_y=new float[descriptor.nnodes];

  rhs=new double[descriptor.nnodes];
  
  combo_gx.allocate(combo.nfields);
  combo_gy.allocate(combo.nfields);
  
  for(k=0;k<combo.nfields;k++) {
    float *topo_x,*topo_y;
/*------------------------------------------------------------------------------
    load data */
    combo.fields[k].acquire(0.);
    size_t count __attribute__((unused))=occurence(combo.fields[k].mask, combo.fields[k].x, combo.fields[k].grid->Hsize());
/*------------------------------------------------------------------------------
    compute gradients */
    combo_gx.fields[k]=combo.fields[k];
    combo_gy.fields[k]=combo.fields[k];
    size=combo.fields[k].nvalues();
    topo_x=new float[size];
    topo_y=new float[size];
    status=map_gradient(*combo.fields[k].grid,combo.fields[k].grid->nx,combo.fields[k].x,combo.fields[k].mask, GEOCENTRIC, topo_x, topo_y);
    combo_gx.fields[k].x=topo_x;
    combo_gy.fields[k].x=topo_y;
/*------------------------------------------------------------------------------
    needed by interpolation routine */
    set_grid_list(combo.fields[k].grid);
    }
  
//  status=map_resolution( topogrid,  &dx,  &dy,  &dz);

  printf("#################################################################\n");
  printf("interpolate depths\n");
  
  RHS_generic(mesh, descriptor, rhs, combo, missing);
  status = LinearSystem_solve(descriptor.massmatrix,rhs);
  for (n=0;n<descriptor.nnodes;n++) {
    depths[n]=factor*rhs[n];
    }

  discretisation=descriptor.type;

  sprintf(comments[0],"%s bathymetry, unit in meters",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, optimal estimate from %s",dataset);
  sprintf(filename,"topo-%s-1.s2r",discretisation_name(discretisation));
  status=quoddy_saver1(filename, descriptor.nnodes, depths,comments);
  
  if(output!=0) {
    for(int n=0;n<mesh.nvtxs;n++) mesh.vertices[n].h=depths[n];
    status=fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID, mesh);
    }

//  rms=CHK_generic(mesh, descriptor, rhs, topogrid, topo, topomask);
//   sprintf(varname,"rmsH-%s",discretisation_name(discretisation));
//   status=archiving_UGdummy2D("mesh-topo.check.nc", mesh, varname, rms,(int) LGP0);

  printf("#################################################################\n");
  printf("interpolate depth's gradients\n");
  RHS_generic(mesh, descriptor, rhs, combo_gx, missing);
  status = LinearSystem_solve(descriptor.massmatrix,rhs);
  for (n=0;n<descriptor.nnodes;n++) {
    slope_x[n]=(float) rhs[n];
    slope_x[n]*=1000.0;
    }
  
  RHS_generic(mesh, descriptor, rhs, combo_gy, missing);
  status = LinearSystem_solve(descriptor.massmatrix,rhs);
  for (n=0;n<descriptor.nnodes;n++) {
    slope_y[n]=(float) rhs[n];
    slope_y[n]*=1000.0;
    }

  sprintf(comments[0],"%s bathymetry slope, dimensionless",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, optimal estimate from %s",dataset);
  sprintf(filename,"slope-%s-1.v2r",discretisation_name(discretisation));
  status=quoddy_saver2(filename, descriptor.nnodes, slope_x, slope_y,comments);

  delete[] rhs;

  delete[] depths;

  delete[] slope_x;
  delete[] slope_y;

  delete[] comments[0];
  delete[] comments[1];

  printf("%s completed (optimal)\n",__func__);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_TOPO(SGfield_t<float> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k, status;
  int nomask=0;
  grid_t *topogrid=new grid_t;
  float *topo,topomask,scale=1.0;
  filestream_t<float> *stream;
  bool debug=false;
  
  stream=(filestream_t<float> *) field.stream;

  printf("#################################################################\n");
  printf("load bathymetry database\n");
//  status=map_loadfield(stream->filename.c_str(), (char *) 0, topogrid, &topo, &topomask);
  status=topo_loadfield(stream->filename.c_str(),stream->varname.c_str(), topogrid, &topo, &topomask, debug);
  
  if(status!=0) return(status);

  if(scale==-1.) {
    printf("#################################################################\n");
    printf("apply -1 factor to bathymetry\n");
    for(j=0;j<topogrid->ny;j++)
      for(i=0;i<topogrid->nx;i++) {
        k=j*topogrid->nx+i;
        if (topo[k]!=topomask) {
/*------------------------------------------------------------------------------
          set ocean positive depth*/
          topo[k]=-topo[k];
          }
        else {
          topo[k]=0;
          }
        }
      }

/*------------------------------------------------------------------------------
  quick fix for limited bathymetry; useful ??? */
  if(nomask==1) {
    printf("#################################################################\n");
    printf("change masked values to zero\n");
    for(j=0;j<topogrid->ny;j++)
      for(i=0;i<topogrid->nx;i++) {
        k=j*topogrid->nx+i;
        if (topo[k]==topomask) {
          topo[k]=0;
          }
        }
      }

  field.x=topo;
  field.mask=topomask;
  field.grid=topogrid;
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int field_interpolation( SGfield_t<float> & field, double x, double y, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0,verbose=0;
  int64_t m;
  
//   status=map_interpolation(*field.grid,field.x,field.mask,x,y,z);
      
  index_interpolation(*field.grid,x,y,&m,field.x,field.mask,z,verbose);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int create_combo(Scombo_t & combo, vector<string> inputs, string name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  int status;
  double t=CNES0jd; /* dummy date for filestream->finalize(...) */
  vector<string> variable;
  
  
  variable.push_back("z");
  variable.push_back("bathymetry");
  variable.push_back("elevation");
  if(name!="") variable.push_back(name);
  
  combo.allocate(inputs.size());
  
  for(n=0;n<combo.nfields;n++) {
    string path="",name=inputs[n];
    
    filestream_t<float> *filestream=new filestream_t<float> (path.c_str(),name.c_str(), 0, 0);
    for(int s=0;s<variable.size();s++) {
      status=filestream->finalize(t,variable[s].c_str());
      if(status==0) break;
      }
    if(status!=0) {
      printf("%s : could not load %s\n", __func__, filestream->filename.c_str());
      for(int s=0;s<variable.size();s++) printf("attemptive variable name %d %s\n",s,variable[s].c_str());
      return(status);
      }
    combo.fields[n].stream=(datastream_t<float> *) filestream;
    combo.fields[n].processor=stream_TOPO;
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int i,j,k,m,n;
  int paire;
  char *keyword,*zone=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*spaire=NULL;
  char *discretisation=NULL;
  mesh_t mesh;
  grid_t topogrid;
  float *topo, topomask, scale=1.0, factor;
  float missing=99999.9;
  int nomask=0;
  discretisation_t descriptorQ0, descriptorQ1;
  vector<string> inputs;
  string variable;
  size_t count;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  bottom topography use:
  ----------------------
  
    Tradionnally models expect positive depths, so in the usual case where bottom
    topography databases provide negative depths, sign must be changed in the
    final depth vector. So it is by deault in mesh-topo
  
    gradient are expexted to be computed on (true) negative depths
    
  other uses:
  -----------

    likely no change in signs will be necessary, so use <--change-sign=no> in
    the line command otions
    
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

bool change_sign=true;
  
  Scombo_t combo;
  bool debug=false;
  
  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"--change-sign=no")==0) {
      change_sign=false;
      n++;
      continue;
      }
    if(strcmp(keyword,"--change-sign=yes")==0) {
      change_sign=true;
      n++;
      continue;
      }
    if(strcmp(keyword,"-nomask")==0) {
      nomask=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-missing")==0) {
      int nitems;
      nitems=sscanf(argv[n+1], "%f",&missing);
      n++;
      n++;
      continue;
      }
    if(strcmp(keyword,"-debug")==0) {
      debug=true;
      n++;
      continue;
      }
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

        case 'd' :
          discretisation= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          if(depthfile==NULL) depthfile = strdup(argv[n+1]);
          inputs.push_back((string) argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          spaire= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
          variable=argv[n+1];
          n++;
          n++;
          break;

        case 'h' :
          print_help(argv[0]);
          n++;
          break;

        default:
          print_help(argv[0]);
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        print_help(argv[0]);
        STDOUT_BASE_LINE("unknown option %s\n",keyword);
        exit(-1);
      }
    free(keyword);
    }

  fe_integrale_init();

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load and initialize mesh structure
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("load mesh and initialize related tables\n");fflush(stdout);
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status!=0) TRAP_ERR_EXIT(status,"unable to read the original mesh in %s (error code %d)\n",meshfile,status);
/**----------------------------------------------------------------------------
  test - test - test - test - test - test - test - test - test - test - test */
//     status=fe_list(&mesh);
    status=quadrangle_list ( mesh, 3, 4);
    if(status!=0) TRAP_ERR_EXIT(status,"unable to build the element list from the original mesh (error code %d)\n",status);
    }
  else {
    STDOUT_BASE_LINE("*** please specify a mesh file with option -m ***\n");
    print_help(argv[0]);
    exit(-1);
    }

  if(spaire==0) {
    STDOUT_BASE_LINE("*** please specify a computational paire with option -p ***\n");
    print_help(argv[0]);
    exit(-1);
    }
  
  if(mesh.ntriangles !=0) {
    status=fe_geometry(&mesh);
    status= fe_edgetable(&mesh,0,0);
    }
  else {
    status= fe_edgetable_Q(&mesh,0);
    for(m=0;m<mesh.nquadrangles;m++) {
      status=fe_initaffine_spherical(mesh, &(mesh.quadrangles[m]), m);
      }
    }
  if(status!=0) TRAP_ERR_EXIT(status,"unable to build the element list from the original mesh (error code %d)\n",status);

  status=create_combo(combo, inputs, variable);
  if(status!=0) TRAP_ERR_EXIT(status,"database error, code %d\n",status);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  load and initialize bathymetry
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("load bathymetry database\n");fflush(stdout);
  status=topo_loadfield(depthfile, variable.c_str(), &topogrid, &topo, &topomask, debug);
  NC_TRAP_ERROR(wexit,status,1,"topo_loadfield() error");
  
  count=occurence(topomask, topo, topogrid.Hsize());

  if(scale==-1.) {
    printf("#################################################################\n");
    printf("apply -1 factor to bathymetry\n");fflush(stdout);
    for(j=0;j<topogrid.ny;j++)
      for(i=0;i<topogrid.nx;i++) {
        k=j*topogrid.nx+i;
        if (topo[k]!=topomask) {
/*------------------------------------------------------------------------------
          set ocean positive depth*/
          topo[k]=-topo[k];
          }
        else {
          topo[k]=0;
          }
        }
      }

/*------------------------------------------------------------------------------
  quick fix for limited bathymetry*/
  if(nomask==1) {
    printf("#################################################################\n");
    printf("change masked values to zero\n");fflush(stdout);
    for(j=0;j<topogrid.ny;j++)
      for(i=0;i<topogrid.nx;i++) {
        k=j*topogrid.nx+i;
        if (topo[k]==topomask) {
          topo[k]=0;
          }
        }
      }


  if(strcmp(spaire,"LGP0xLGP1")==0) {
    paire=LGP0xLGP1;
    }
  else if(strcmp(spaire,"DGP1xLGP2")==0) {
    paire=DGP1xLGP2;
    }
  else if(strcmp(spaire,"DNP1xLGP2")==0) {
    paire=DNP1xLGP2;
    }
  else if(strcmp(spaire,"Q1xQ0")==0) {
    paire=CQP1xCQP0;
    }
  else if(strcmp(spaire,"CQN1xCQQ0")==0) {
    paire=CQP1xCQP0;
    }
  else {
    STDOUT_BASE_LINE("*** computational paire not recognised : give one of LGP0xLGP1, DGP1xLGP2, DNP1xLGP2, Q1xQ0 or CQN1xCQQ0 ***\n");
    exit(-1);
    }
  
  
  printf("#################################################################\n");
  printf("derive model bathymetry\n");fflush(stdout);
  switch(paire) {
    case LGP0xLGP1:
      topogrid.free();
      delete[] topo;
      status= discretisation_init(&mesh, LGP0);
      status= discretisation_init(&mesh, LGP1);
      if(change_sign) factor=-1.0;
      else factor=1.;
      missing*=factor;
      status= fe_topo01_generic(mesh, mesh.LGP1descriptor, combo, depthfile, factor, missing);
      status= fe_topo02_generic(mesh, mesh.LGP1descriptor, combo, depthfile, factor, missing, output);
      break;
    case DGP1xLGP2:
      printf("Warning : sign change option not managed\n");
      status= discretisation_init(&mesh, LGP0);
      status= fe_topo01_generic_obsolete(mesh, mesh.LGP0descriptor, topogrid, topo, topomask, depthfile);
      status= fe_topo02_generic_obsolete(mesh, mesh.LGP0descriptor, topogrid, topo, topomask, depthfile);

      status= discretisation_init(&mesh, DGP2);
      status= fe_topo01_generic_obsolete(mesh, mesh.DGP2descriptor, topogrid, topo, topomask, depthfile);
      status= fe_topo02_generic_obsolete(mesh, mesh.DGP2descriptor, topogrid, topo, topomask, depthfile);

      status= discretisation_init(&mesh, LGP1);
      status= fe_topo01_generic_obsolete(mesh, mesh.LGP1descriptor, topogrid, topo, topomask, depthfile);
      status= fe_topo02_generic_obsolete(mesh, mesh.LGP1descriptor, topogrid, topo, topomask, depthfile);

      status= discretisation_init(&mesh, LGP2);
      status= fe_topo01_generic_obsolete(mesh, mesh.LGP2descriptor, topogrid, topo, topomask, depthfile);
      status= fe_topo02_generic_obsolete(mesh, mesh.LGP2descriptor, topogrid, topo, topomask, depthfile);

      status= discretisation_init(&mesh, DGP1);
      status= fe_topo01_generic_obsolete(mesh, mesh.DGP1descriptor, topogrid, topo, topomask, depthfile);
      status= fe_topo02_generic_obsolete(mesh, mesh.DGP1descriptor, topogrid, topo, topomask, depthfile);
      break;
    case DNP1xLGP2:
      printf("Warning : sign change option not managed\n");
      status= discretisation_init(&mesh, DNP1);
      status= fe_topo01_generic_obsolete(mesh, mesh.DNP1descriptor, topogrid, topo, topomask, depthfile);
      status= fe_topo02_generic_obsolete(mesh, mesh.DNP1descriptor, topogrid, topo, topomask, depthfile);
      status= discretisation_init(&mesh, LGP2);
      status= fe_topo01_generic_obsolete(mesh, mesh.LGP2descriptor, topogrid, topo, topomask, depthfile);
      status= fe_topo02_generic_obsolete(mesh, mesh.LGP2descriptor, topogrid, topo, topomask, depthfile);

//       status= discretisation_init(&mesh, DNP1);
//       status= fe_topo01_generic_obsolete(mesh, mesh.DNP1descriptor, topogrid, topo, topomask, depthfile);
//       status= fe_topo02_generic(mesh, mesh.DNP1descriptor, topogrid, topo, topomask, depthfile);
      break;
    case CQP1xCQP0:
      printf("Warning : sign change option not managed\n");
      status= discretisation_init(&mesh, CQP0);
      descriptorQ1=initCQP1(mesh, 0);
      descriptorQ1.massmatrix.ordering=new ordering_t();
      status=dMassMatrix(mesh, descriptorQ1, &descriptorQ1.massmatrix);
      status= fe_topo01_generic_obsolete(mesh, descriptorQ1, topogrid, topo, topomask, depthfile);
      status= fe_topo02_generic_obsolete(mesh, descriptorQ1, topogrid, topo, topomask, depthfile);
      descriptorQ0=initCQP0(mesh, 0);
      descriptorQ0.massmatrix.ordering=new ordering_t();
      status=dMassMatrix(mesh, descriptorQ0, &descriptorQ0.massmatrix);
      status= fe_topo01_generic_obsolete(mesh, descriptorQ0, topogrid, topo, topomask, depthfile);
      status= fe_topo02_generic_obsolete(mesh, descriptorQ0, topogrid, topo, topomask, depthfile);
      break;
    default:
      break;
    }
  
  STDOUT_BASE_LINE("end of mesh-topo ... \n");
  exit(0);
}
