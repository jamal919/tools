
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

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
#include "discretisation.h"
#include "mass.h"
#include "matrix.h"
#include "netcdf-classes.h"
#include "netcdf-proto.h"
#include "datastream.h"

#include "map.def"

#include "zapper.h"

extern void fe_LGP0_to_LGP1(double *in, double *out, mesh_t mesh);

extern  discretisation_t initCQP1(mesh_t mesh, int type);
extern  discretisation_t initCQP0(mesh_t mesh, int type);

int field_interpolation( Sfield_t<float> & field, double x, double y, float *z);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_edgetable_quick(mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,j1,j2,k,l,m,n1,n2,n;
  int   nedges,chk=0,count=0,tmp=0;
  int   nelt,nndes;
  int   interior=0,boundary=0,weird=0;
  int   ring[4]={0,1,2,0};
  int   *connexions;
  double t1,t2,p1,p2;
  edge_t *edges;
  triangle_t *elt,*ptr;

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
        if(mesh->type==0) t2=geo_recale(t2,t1,(double) 180.0);
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
  int   k,l,m,i,j,j1,j2,n1,n2,n,status;
  int   nedges,chk=0,count=0,tmp=0;
  int   nelts,nndes;
  int   interior=0,boundary=0,weird=0;
  int   ring[4]={0,1,2,0};
  int   *ncells,ncellsmax=0,**cells;
  edge_t *edges;
  triangle_t *elt,*ptr;
  double t1,t2,p1,p2,h1,h2;
  mesh_t work;

  elt=mesh.triangles;

  nelts =mesh.ntriangles;
  nndes =mesh.nvtxs;
  nedges=mesh.nedges;
  edges =mesh.edges;

  status=fe_init(&work);
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
//     if(mesh.type==0) t2=geo_recale(t2,t1,(double) 180.0);
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
//    status=fe_initaffine(&work,work.ntriangles);
    work.ntriangles++;
    work.triangles[work.ntriangles].vertex[0]=elt[k].vertex[1];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[0];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[2];
//    status=fe_initaffine(&work,work.ntriangles);
    work.ntriangles++;
    work.triangles[work.ntriangles].vertex[0]=elt[k].vertex[2];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[1];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[0];
//    status=fe_initaffine(&work,work.ntriangles);
    work.ntriangles++;
    work.triangles[work.ntriangles].vertex[0]=nndes+elt[k].edges[2];
    work.triangles[work.ntriangles].vertex[1]=nndes+elt[k].edges[0];
    work.triangles[work.ntriangles].vertex[2]=nndes+elt[k].edges[1];
//    status=fe_initaffine(&work,work.ntriangles);
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
  int   i,j,k,l,m,n,status;
  int   count=0,n1,n2;
  double t1,p1,t2,p2,lon,lat;
  int   *used;
  triangle_t *e;
  vertex_t  *v;
  mesh_t work;

  status=fe_init(&work);

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
    t2=geo_recale(t2,t1,(double) 180.0);
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

  status=fe_e2n(&work);

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
  int   i,j,k,l,m,n,status;
  int   count,kmax=3;
  double t,p,rr,rate,w;
  float h;
  mesh_t work,refined[10];

  rr=0;
  for (i=0;i<3;i++) {
    n=mesh.triangles[target].edges[i];
    rr=MAX(rr,mesh.edges[n].L);
    }

  rate=rr/r;
  kmax=0;
  while(rate > 0.5) {
    rate=rate/2.0;
    kmax++;
    }

  refined[0]=fe_splitelement(mesh, target);
  status=fe_edgetable(&(refined[0]),0,0);
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
  int   i,j,k,l,m,n,status;
  int   count,kmax=3;
  double t,p,rr,rate,w;
  float h;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;

  gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);

  rr=0;
  for (i=0;i<3;i++) {
    n=mesh.triangles[target].edges[i];
    rr=MAX(rr,mesh.edges[n].L);
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
        printf("trouble\n");
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
  int count, i, j, k, l, m, n, status;
  int n1, n2;
  int ni, nj, row, column;
  double C, area;
  double beta[3];
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double sum, **base;
  float h;
  int discretisation=descriptor.id;
  
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
          printf("trouble\n");
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
  int count, i, j, k, l, m, n, status;
  int n1, n2;
  int ni, nj, row, column;
  double C, area;
  double beta[3];
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double sum, **base;
  float h;
  int discretisation=descriptor.id;
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
          printf("trouble\n");
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

  void RHS_generic_T(mesh_t mesh, discretisation_t descriptor, double *rhs, Scombo_t combo)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Triangle Gauss integration

----------------------------------------------------------------------*/
  int count, i, j, k, l, m, n, status;
  int n1, n2;
  int ni, nj, row, column;
  double C, area;
  double beta[3];
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double sum, **base;
  float h,hh;
  int discretisation=descriptor.id;
  float mask;
  
  status=gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);
  
  for(n = 0; n < descriptor.nnodes; n++)
    rhs[n] = 0;
  
  base=new double*[gauss_n];
  for (k=0;k<gauss_n;k++) {
    base[k]=new double[descriptor.nnpe];
    status=fe_LGPbase(mesh, discretisation, gauss_x[k], gauss_y[k], base[k]);
    }

  mask=combo.fields[0].mask;
  
  for(m = 0; m < mesh.ntriangles; m++) {
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
          if(hh!=combo.fields[l].mask) {
            h=hh;
            }
//           else {
//             printf("trouble\n");
//             }
          }
        if(h==mask) {
          printf("trouble\n");
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

  void RHS_generic_Q(mesh_t mesh, discretisation_t descriptor, double *rhs, Scombo_t combo)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *---------------------------------------------------------------------------

  Quadrangle Gauss integration

-----------------------------------------------------------------------------*/
  int count, i, j, k, l, m, n, status;
  int n1, n2;
  int ni, nj, row, column;
  double C, area;
  double beta[3];
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double sum, **base;
  float h;
  int discretisation=descriptor.id;
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
            printf("trouble\n");
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

  void RHS_generic(mesh_t mesh, discretisation_t descriptor, double *rhs, Scombo_t combo)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
  switch(mesh.ntriangles) {
    case 0:
      RHS_generic_Q(mesh, descriptor, rhs, combo);
      break;
    default:
      RHS_generic_T(mesh, descriptor, rhs, combo);
      break;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *CHK_generic(mesh_t mesh, discretisation_t descriptor, double *UGtopo, grid_t topogrid, float *buffers,float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Mass matrix for P2 integrale scalar product

----------------------------------------------------------------------*/
  int count, i, j, k, l, m, n, status;
  int n1, n2;
  int ni, nj, row, column;
  double C, area;
  double beta[3];
  double t,p;
  double gauss_x[16], gauss_y[16], gauss_w[16];
  int    gauss_n=16;
  double sum, **base;
  gauss_init(gauss_n,gauss_x,gauss_y,gauss_w);
  float h,*UGh;
  double *rms;
  double global;
  int discretisation=descriptor.id;
    
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
        printf("trouble\n");
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

  int fe_topo01_generic_obsolete(mesh_t mesh, discretisation_t descriptor, grid_t topogrid, float *topo, float topomask, char *dataset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,k,l,m,n,status;
  double t,p;
  float h;
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

  for (n=0;n<descriptor.nnodes;n++) {
    t=descriptor.nodes[n].lon;
    p=descriptor.nodes[n].lat;
    t=map_recale(topogrid,t);
    status=map_interpolation(topogrid,topo,topomask,t,p,&(depths[n]));
    if(depths[n]==topomask) {
      printf("trouble\n");
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
    
  discretisation=descriptor.id;
  
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

  printf("fe_topo01_generic_obsolete completed (raw interpolation)\n\n");

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_topo02_generic_obsolete(mesh_t mesh, discretisation_t descriptor, grid_t topogrid, float *topo,float topomask, char *dataset)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,k,l,m,n,status;
  int   nbuffers;
  double t,p,tref;
  double *M;
  ordering_t *ordering;
  float h,*depths,*slope_x,*slope_y;
  double **buffer_LGP0,*buffer_LGP1,*rhs;
  float **buffers;
  float *topo_x,*topo_y;
  size_t size;
  char *comments[2], *filename, *varname;
  mesh_t work;
  double *mean;
  double *dx,*dy,*dz,r;
  int discretisation;
  double *rms;

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

  discretisation=descriptor.id;

  sprintf(comments[0],"%s bathymetry, unit in meters",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, optimal estimate from %s",dataset);
  sprintf(filename,"topo-%s-1-old.s2r",discretisation_name(discretisation));
  status=quoddy_saver1(filename, descriptor.nnodes, depths,comments);
  
  rms=CHK_generic(mesh, descriptor, rhs, topogrid, topo, topomask);
  sprintf(varname,"rmsH-%s",discretisation_name(discretisation));
  status=archiving_UGdummy2D("mesh-topo.check.nc", mesh, varname, rms,(int) LGP0);

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

  printf("fe_topo02_generic_obsolete completed (optimal)\n\n");
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_topo01_generic(mesh_t mesh, discretisation_t descriptor, Scombo_t combo, char *dataset, float factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,k,l,m,n,status;
  double t,p;
  float h;
  float *depths,*slope_x,*slope_y,mask=99999.9,z;
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
  for (n=0;n<descriptor.nnodes;n++) depths[n]=mask;
  
  for(k=0;k<combo.nfields;k++) {
    combo.fields[k].acquire(0.);
    for (n=0;n<descriptor.nnodes;n++) {
      t=descriptor.nodes[n].lon;
      p=descriptor.nodes[n].lat;
      t=map_recale(*combo.fields[k].grid,t);
//      status=map_interpolation(*combo.fields[k].grid,combo.fields[k].grid->nx,combo.fields[k].x,combo.fields[k].mask,t,p,&z);
      status=field_interpolation(combo.fields[k],t,p,&z);
      if(z==combo.fields[k].mask) {
//        printf("trouble\n");
        }
/**----------------------------------------------------------------------------
      model needs positive depths */
      if(z!=combo.fields[k].mask) {
        depths[n]=factor*z;
        }
      }
    }
    
  for (n=0;n<descriptor.nnodes;n++) {
    if(depths[n]==mask) {
       printf("trouble\n");
       }
    }
   
  discretisation=descriptor.id;
  
  sprintf(comments[0],"%s bathymetry, unit in meters",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, direct interpolation from %s",dataset);
  sprintf(filename,"topo-%s-0.s2r",discretisation_name(discretisation));
  status=quoddy_saver1(filename, descriptor.nnodes, depths,comments);

//   rms=CHK_generic(mesh, descriptor, rhs, topogrid, topo, topomask);
//   sprintf(varname,"rmsH-%s",discretisation_name(discretisation));
//   status=archiving_UGdummy2D("mesh-topo.check.nc", mesh, varname, rms,(int) LGP0);


  printf("#################################################################\n");
  printf("interpolate depth's gradients\n");
  
  slope_x=new float[descriptor.nnodes];
  slope_y=new float[descriptor.nnodes];
  
  for(k=0;k<combo.nfields;k++) {
    size=combo.fields[k].nvalues();
    topo_x=new float[size];
    topo_y=new float[size];
    status= map_gradient(*combo.fields[k].grid,combo.fields[k].grid->nx,combo.fields[k].x,combo.fields[k].mask, GEOCENTRIC, topo_x, topo_y);

    for (n=0;n<descriptor.nnodes;n++) {
      t=descriptor.nodes[n].lon;
      p=descriptor.nodes[n].lat;
      t=map_recale(*combo.fields[k].grid,t);
      status=map_interpolation(*combo.fields[k].grid,combo.fields[k].grid->nx,topo_x,combo.fields[k].mask,t,p,&z);
//       if(z==combo.fields[k].mask) {
//         printf("trouble\n");
//         }
      if(z!=combo.fields[k].mask) {
        slope_x[n]=z*1000.0;
        }
      status=map_interpolation(*combo.fields[k].grid,combo.fields[k].grid->nx,topo_y,combo.fields[k].mask,t,p,&z);
//       if(z==combo.fields[k].mask) {
//         printf("trouble\n");
//         }
      if(z!=combo.fields[k].mask) {
        slope_y[n]=z*1000.0;
        }
      }
      
    delete[] topo_x;
    delete[] topo_y;
    }

  for (n=0;n<descriptor.nnodes;n++) {
    if(slope_x[n]==mask) {
       printf("trouble\n");
       }
    if(slope_y[n]==mask) {
       printf("trouble\n");
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

  printf("fe_topo01_generic completed (raw interpolation)\n");

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_topo02_generic(mesh_t mesh, discretisation_t descriptor, Scombo_t combo, char *dataset, float factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,k,l,m,n,status;
  int   nbuffers;
  double t,p,tref;
  double *M;
  ordering_t *ordering;
  float h,*depths,*slope_x,*slope_y;
  double **buffer_LGP0,*buffer_LGP1,*rhs;
  float **buffers;
//  float *topo_x,*topo_y;
  size_t size;
  char *comments[2], *filename, *varname;
  mesh_t work;
  double *mean;
  double *dx,*dy,*dz,r;
  int discretisation;
  double *rms;
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
    combo.fields[k].acquire(0.);
    combo_gx.fields[k]=combo.fields[k];
    combo_gy.fields[k]=combo.fields[k];
    size=combo.fields[k].nvalues();
    topo_x=new float[size];
    topo_y=new float[size];
    status= map_gradient(*combo.fields[k].grid,combo.fields[k].grid->nx,combo.fields[k].x,combo.fields[k].mask, GEOCENTRIC, topo_x, topo_y);
    combo_gx.fields[k].x=topo_x;
    combo_gy.fields[k].x=topo_y;
    }
    
//  status=map_resolution( topogrid,  &dx,  &dy,  &dz);

  printf("#################################################################\n");
  printf("interpolate depths\n");
  
  RHS_generic(mesh, descriptor, rhs, combo);
  status = LinearSystem_solve(descriptor.massmatrix,rhs);
  for (n=0;n<descriptor.nnodes;n++) {
    depths[n]=factor*rhs[n];
    }

  discretisation=descriptor.id;

  sprintf(comments[0],"%s bathymetry, unit in meters",discretisation_name(discretisation));
  sprintf(comments[1],"Created by mesh-topo, optimal estimate from %s",dataset);
  sprintf(filename,"topo-%s-1.s2r",discretisation_name(discretisation));
  status=quoddy_saver1(filename, descriptor.nnodes, depths,comments);
  
//  rms=CHK_generic(mesh, descriptor, rhs, topogrid, topo, topomask);
//   sprintf(varname,"rmsH-%s",discretisation_name(discretisation));
//   status=archiving_UGdummy2D("mesh-topo.check.nc", mesh, varname, rms,(int) LGP0);

  printf("#################################################################\n");
  printf("interpolate depth's gradients\n");
  RHS_generic(mesh, descriptor, rhs, combo_gx);
  status = LinearSystem_solve(descriptor.massmatrix,rhs);
  for (n=0;n<descriptor.nnodes;n++) {
    slope_x[n]=(float) rhs[n];
    slope_x[n]*=1000.0;
    }
    
  RHS_generic(mesh, descriptor, rhs, combo_gy);
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

//   delete[] topo_x;
//   delete[] topo_y;

  delete[] slope_x;
  delete[] slope_y;

  delete[] comments[0];
  delete[] comments[1];

  printf("fe_topo02_generic completed (optimal)\n");
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_TOPO(Sfield_t<float> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n, status;
  int nomask=0;
  grid_t *topogrid=new grid_t;
  float *topo,topomask,h,scale=1.0;
  filestream_t<float> *stream;
  
  stream=(filestream_t<float> *) field.stream;

  printf("#################################################################\n");
  printf("load bathymetry database\n");
  status=map_loadfield(stream->filename.c_str(), (char *) 0, topogrid, &topo, &topomask);

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
   quick fix for limited bathymetry*/
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

int field_interpolation( Sfield_t<float> & field, double x, double y, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=map_interpolation(*field.grid,field.grid->nx,field.x,field.mask,x,y,z);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int create_combo(Scombo_t & combo, vector<string> inputs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  int status;
//  Scombo combo;
  double t;
  
  combo.allocate(inputs.size());
  
  for(n=0;n<combo.nfields;n++) {
    string path="",name=inputs[n];
    
    filestream_t<float> *filestream=new filestream_t<float> (path.c_str(),name.c_str());
    status=filestream->finalize(t,"z");
    
    combo.fields[n].stream=(datastream_t<float> *) filestream;
    combo.fields[n].processor=stream_TOPO;
    }
      
//       combo.fields[0].descriptor=get_descriptor_address(mesh, u2D_discretisation);
//       combo.fields[1].descriptor=get_descriptor_address(mesh, u2D_discretisation);
//       combo.fields[0].allocate();
//       combo.fields[1].allocate();
//       combo.stream=Ustream;
//       combo.processor=polar_conversion;
//       }
//     target=t;
//     for(int k=0;k<atmosphere.winds.nframes;k++) {
//       combo.acquire(target);
//       target=combo.next;
//       }
//     }
 
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option,paire;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*spaire=NULL;
  char *discretisation=NULL;
  mesh_t mesh;
  int *selected,*targeted;
  grid_t topogrid;
  float *topo,topomask,h,scale=1.0;
  int nomask=0;
  discretisation_t descriptorQ0, descriptorQ1;
  vector<string> inputs;
  Scombo_t combo;
  int nfiles;
  char **filenames;
  
  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    if(strcmp(keyword,"-inv")==0) {
      scale=-1;
      n++;
      continue;
      }
    if(strcmp(keyword,"-nomask")==0) {
      nomask=1;
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

  fe_integrale_init();

/* *----------------------------------------------------------------------
  load and initialize mesh structure*/
  printf("#################################################################\n");
  printf("load mesh ant initialize related tables\n");
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status!=0) {
      printf("unable to read the original mesh in %s\n",meshfile);
      goto error;
      }
/**----------------------------------------------------------------------------
  test - test - test - test - test - test - test - test - test - test - test */
    status=fe_list(&mesh);
//    status=quadrangle_list ( mesh);
    if(status!=0) {
      printf("unable to build the element list from the original mesh\n");
      goto error;
      }
    }
 else {
   printf("no mesh file specified; abort...\n");
   goto error;
   }

  if(mesh.ntriangles !=0) {
    status=fe_edgetable(&mesh,0,0);
    }
  else {
    status=fe_edgetable_Q(&mesh,0);
    for(m=0;m<mesh.nquadrangles;m++) {
      status=fe_initaffine_spherical(&mesh, &(mesh.quadrangles[m]), m);
      }
    }
  if(status!=0) {
    printf("unable to build the edge list from the original mesh\n");
    goto error;
    }

  nfiles=inputs.size();
//   filenames=new char*[nfiles];
//   filenames[0]=strdup(depthfile);
//  status=create_combo(combo, nfiles, filenames);
  status=create_combo(combo, inputs);
  
  printf("#################################################################\n");
  printf("load bathymetry database\n");
  status=map_loadfield(depthfile, (char *) 0,&topogrid,&topo,&topomask);

  if(scale==-1.) {
    printf("#################################################################\n");
    printf("apply -1 factor to bathymetry\n");
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
    printf("change masked values to zero\n");
    for(j=0;j<topogrid.ny;j++)
      for(i=0;i<topogrid.nx;i++) {
        k=j*topogrid.nx+i;
        if (topo[k]==topomask) {
          topo[k]=0;
          }
        }
      }


  if(spaire==0) {
    printf("computational paire not specified, abort...\n");
    goto error;
    }
  else if(strcmp(spaire,"LGP0xLGP1")==0) {
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
  else {
    goto error;
    }
  
  
  printf("#################################################################\n");
  printf("derive model bathymetry\n");
  switch(paire) {
    case LGP0xLGP1:
      status= discretisation_init(&mesh, LGP0);
      status= discretisation_init(&mesh, LGP1);
//       status=fe_topo01_generic_obsolete(mesh, mesh.LGP1descriptor, topogrid, topo, topomask, depthfile);
//       status=fe_topo02_generic_obsolete(mesh, mesh.LGP1descriptor, topogrid, topo, topomask, depthfile);
//       delete[] topo;
      status=fe_topo01_generic(mesh, mesh.LGP1descriptor, combo, depthfile);
      status=fe_topo02_generic(mesh, mesh.LGP1descriptor, combo, depthfile);
      break;
    case DGP1xLGP2:
      status= discretisation_init(&mesh, LGP0);
      status=fe_topo01_generic_obsolete(mesh, mesh.LGP0descriptor, topogrid, topo, topomask, depthfile);
      status=fe_topo02_generic_obsolete(mesh, mesh.LGP0descriptor, topogrid, topo, topomask, depthfile);

      status= discretisation_init(&mesh, DGP2);
      status=fe_topo01_generic_obsolete(mesh, mesh.DGP2descriptor, topogrid, topo, topomask, depthfile);
      status=fe_topo02_generic_obsolete(mesh, mesh.DGP2descriptor, topogrid, topo, topomask, depthfile);

      status= discretisation_init(&mesh, LGP1);
      status=fe_topo01_generic_obsolete(mesh, mesh.LGP1descriptor, topogrid, topo, topomask, depthfile);
      status=fe_topo02_generic_obsolete(mesh, mesh.LGP1descriptor, topogrid, topo, topomask, depthfile);

      status= discretisation_init(&mesh, LGP2);
      status=fe_topo01_generic_obsolete(mesh, mesh.LGP2descriptor, topogrid, topo, topomask, depthfile);
      status=fe_topo02_generic_obsolete(mesh, mesh.LGP2descriptor, topogrid, topo, topomask, depthfile);

      status= discretisation_init(&mesh, DGP1);
      status=fe_topo01_generic_obsolete(mesh, mesh.DGP1descriptor, topogrid, topo, topomask, depthfile);
      status=fe_topo02_generic_obsolete(mesh, mesh.DGP1descriptor, topogrid, topo, topomask, depthfile);
      break;
    case DNP1xLGP2:
      status= discretisation_init(&mesh, DNP1);
      status=fe_topo01_generic_obsolete(mesh, mesh.DNP1descriptor, topogrid, topo, topomask, depthfile);
      status=fe_topo02_generic_obsolete(mesh, mesh.DNP1descriptor, topogrid, topo, topomask, depthfile);
      status= discretisation_init(&mesh, LGP2);
      status=fe_topo01_generic_obsolete(mesh, mesh.LGP2descriptor, topogrid, topo, topomask, depthfile);
      status=fe_topo02_generic_obsolete(mesh, mesh.LGP2descriptor, topogrid, topo, topomask, depthfile);

//       status= discretisation_init(&mesh, DNP1);
//       status=fe_topo01_generic_obsolete(mesh, mesh.DNP1descriptor, topogrid, topo, topomask, depthfile);
//       status=fe_topo02_generic(mesh, mesh.DNP1descriptor, topogrid, topo, topomask, depthfile);
      break;
    case CQP1xCQP0:
//      status= discretisation_init(&mesh, DNP1);
      descriptorQ1=initCQP1(mesh, 0);
      descriptorQ1.massmatrix.ordering=new ordering_t();
      status=dMassMatrix(mesh, descriptorQ1, &descriptorQ1.massmatrix);
      status=fe_topo01_generic_obsolete(mesh, descriptorQ1, topogrid, topo, topomask, depthfile);
      status=fe_topo02_generic_obsolete(mesh, descriptorQ1, topogrid, topo, topomask, depthfile);
      descriptorQ0=initCQP0(mesh, 0);
      descriptorQ0.massmatrix.ordering=new ordering_t();
      status=dMassMatrix(mesh, descriptorQ0, &descriptorQ0.massmatrix);
      status=fe_topo01_generic_obsolete(mesh, descriptorQ0, topogrid, topo, topomask, depthfile);
      status=fe_topo02_generic_obsolete(mesh, descriptorQ0, topogrid, topo, topomask, depthfile);
      break;
    default:
      break;
    }

end: __OUT_BASE_LINE__("end of mesh-topo ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
