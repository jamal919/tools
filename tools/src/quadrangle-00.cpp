
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.def"
#include "fe.h"

#include "map.h"
#include "geo.h"
#include "swap.h"
#include "maths.h"
#include "matrix.h"

#include "polygons.hpp"
#include "exceptions.hpp"

#include "polygones.h"

extern  int fe_savemeshNC2D_new(const char *filename, mesh_t & mesh, int option);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_affine_directe(quadrangle_t e, double t, double p,double *ksi,double *eta,int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double teta,phi,dt,dp,tref,x,y;
  double c,u,v;

/*------------------------------------------------------------------------------
  comput ksi,eta coordinates in reference square

  first compute x,y coordinates in reference triangle:

     longitude=t(1)+x*[t(2)-t(1)]+y*[t(3)-t(1)]
     latitude =p(1)+x*[p(2)-p(1)]+y*[p(3)-p(1)]

     determinant = [t(2)-t(1)]*[p(3)-p(1)]-[t(3)-t(1)]*[p(2)-p(1)]

     x = ([t-t(1)]*[p(3)-p(1)]-[t(3)-t(1)]*[p-p(1)])/determinant
       = ([t-t(1)]*cpy-cty*[p-p(1)])/determinant
       
     y = (ctx*[p-p(1)]-[t-t(1)]*cpx)/determinant

     Units : degrees

  then compute u,v coordinates in reference square:

    u= m0 * x / (n0 + n1 * x + n2 * y)
    v= m0 * y / (n0 + n1 * x + n2 * y)

-----------------------------------------------------------------------*/

  tref=e.t_base;

  if(mode==0) teta=geo_recale(t,tref,(double)180.0);
  else teta=t;
  phi = p;

  dt=(teta-e.t_base)*d2r;
  dp=(phi -e.p_base)*d2r;

  x=e.dxdt*dt+e.dxdp*dp;
  y=e.dydt*dt+e.dydp*dp;

  c=(e.n0+e.n1*x+e.n2*y);

  u=e.m0*x/c;
  v=e.m1*y/c;

  *ksi= u;
  *eta= v;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_affine_inverse(mesh_t & mesh, quadrangle_t e, double *t, double *p,double x,double y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double beta[4];
  int i, n;

/*------------------------------------------------------------------------------
  comput ksi,eta coordinates in reference square

  first compute x,y coordinates in reference triangle:

     longitude=t(1)+x*[t(2)-t(1)]+y*[t(3)-t(1)]
     latitude =p(1)+x*[p(2)-p(1)]+y*[p(3)-p(1)]

     determinant = [t(2)-t(1)]*[p(3)-p(1)]-[t(3)-t(1)]*[p(2)-p(1)]

     x = ([t-t(1)]*[p(3)-p(1)]-[t(3)-t(1)]*[p-p(1)])/determinant
       = ([t-t(1)]*cpy-cty*[p-p(1)])/determinant
     y = (ctx*[p-p(1)]-[t-t(1)]*cpx)/determinant

     Units : degrees

  then compute u,v coordinates in reference square:

    u= m0 * x / (n0 + n1 * x + n2 * y)
    v= m0 * y / (n0 + n1 * x + n2 * y)
    
  inverse:
    
    

-----------------------------------------------------------------------*/
  fe_QP1base(x, y, beta);

  *t = 0;
  *p = 0;
  for(i = 0; i < 4; i++) {
    n=e.vertex[i];
    *t += beta[i] * mesh.vertices[n].lon;
    *p += beta[i] * mesh.vertices[n].lat;
    }

  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_initaffine_spherical(const mesh_t & mesh, quadrangle_t *q, int m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------------

  Compute transformation coefficients

  t = s1 t1 +... s4 t4

  p = s1 p1 +... s4 p4


  s1=(1-x)*(1-y)
  s2=   x *(1-y)
  s3=   x * y
  s4=(1-x)* y

  dt/dx =  (1-y) (t2-t1) + y (t3-t4)
  dt/dy =  (1-x) (t4-t1) + x (t3-t2)
  
  !!! not sure about matrix ordering

       | (1-y) (t2-t1) + y (t3-t4)    (1-y) (p2-p1) + y (p3-p4) |
  J =  |                                                        |
       | (1-x) (t4-t1) + x (t3-t2)    (1-x) (p4-p1) + x (p3-p2) |


       | (t2-t1) + y (t1-t2+t3-t4)    (p2-p1) + y (p1-p2+p3-p4) |
  J =  |                                                        |
       | (t4-t1) + x (t1-t4+t3-t2)    (p4-p1) + x (p1-p4+p3-p2) |


       | (t2-t1) + y (t1-t2+t3-t4)    (p2-p1) + y (p1-p2+p3-p4) |
  J =  |                                                        |
       | (t4-t1) + x (t1-t2+t3-t4)    (p4-p1) + x (p1-p2+p3-p4) |


       | a1 + y a2    c1 + y c2 |
  J =  |                                                        |
       | b1 + x b2    d1 + x d2 |


  jacobian = det(J) = (a1 d1 - b1 c1) + x (a1 d2 -c1 b2) + y (a2 d1 - b1 c2) + xy (a2 d2 - b2 c2)
  
  (0,0) : a1 d1 - b1 c1
  
  (1,0) : (a1 d1 - b1 c1) + (a1 d2 -c1 b2)
  
  (0,1) : (a1 d1 - b1 c1) + (a2 d1 - b1 c2)
  
  (1,1) : (a1 d1 - b1 c1) + (a1 d2 -c1 b2) + (a2 d1 - b1 c2) + (a2 d2 - b2 c2)
  
  Square special case :
  
  a1=1 a2=0
  b1=0 b2=0
  c1=0 c2=0
  d1=1 d2=0
  
  
  
  integrale(jacobian) = (a1 d1 - b1 c1) + 0.5 (a1 d2 -c1 b2 + a2 d1 - b1 c2) + 0.25 (a2 d2 - b2 c2)

  http://en.wikipedia.org/wiki/Quadrilateral

  https://www.geometrictools.com/Documentation/PerspectiveMappings.pdf
  

------------------------------------------------------------------------------*/
{
  int n1,n2,n3,n4,status,dum;
  double t1,t2,t3,t4,p1,p2,p3,p4;
  double a1,a2,b1,b2,c1,c2,d1,d2;
/*------------------------------------------------------------------------------
  local earth radius in meters*/
  double R=6.3675e+06;

/*------------------------------------------------------------------------------
  local*/
  double tref,pref;
  double ctx,cty,cpx,cpy,jacobien;
/*------------------------------------------------------------------------------
  passed through triangle structure*/
  double dxdt,dydt,dxdp,dydp;
  vector2D_t u,v,w;
  double c,d;

 redo:

/* *----------------------------------------------------------------------------
  lower-left triangle */
  n1=q->vertex[0];
  n2=q->vertex[1];
  n3=q->vertex[3];
  n4=q->vertex[2];  /// <------- HERE notation not consistent with comments 3<->4

  t1=mesh.vertices[n1].lon;
  t2=mesh.vertices[n2].lon;
  t3=mesh.vertices[n3].lon;
  t4=mesh.vertices[n4].lon;

  if(mesh.type==0) {
    t2=geo_recale(t2,t1,(double) 180.0);
    t3=geo_recale(t3,t1,(double) 180.0);
    t4=geo_recale(t4,t1,(double) 180.0);
    }

  p1=mesh.vertices[n1].lat;
  p2=mesh.vertices[n2].lat;
  p3=mesh.vertices[n3].lat;
  p4=mesh.vertices[n4].lat;

  tref=  t1;
  pref=  p1;

  ctx=  t2-tref;
  cty=  t3-tref;
  cpx=  p2-pref;
  cpy=  p3-pref;

  jacobien=ctx*cpy-cty*cpx;

  if (jacobien > 0.0) {
    dxdt= cpy/(jacobien*d2r);
    dydt=-cpx/(jacobien*d2r);
    dxdp=-cty/(jacobien*d2r);
    dydp= ctx/(jacobien*d2r);
    status=MESH_STATUS_OK;
    }
  else if(jacobien == 0.0) {
    dxdt= 0;
    dydt= 0;
    dxdp= 0;
    dydp= 0;
    status=MESH_STATUS_FLAT_ELEMENT;
    }
  else {
    dum=q->vertex[3];
    q->vertex[3]=q->vertex[1];
    q->vertex[1]=dum;
    status=MESH_STATUS_CW_ELEMENT;
    goto redo;
    }

  q->dxdt  =  dxdt;
  q->dydt  =  dydt;
  q->dxdp  =  dxdp;
  q->dydp  =  dydp;

  q->dtdx  =  ctx;
  q->dpdx  =  cpx;
  q->dtdy  =  cty;
  q->dpdy  =  cpy;

/* *----------------------------------------------------------------------------
  quadrangle */
  u=vector2D_t(t1,p1,t2,p2);
  v=vector2D_t(t1,p1,t3,p3);
  w=vector2D_t(t1,p1,t4,p4);

  w=math_vector_coordinates02(w, u, v);

  c=w.x;
  d=w.y;

  q->m0=d*(1.-c-d);
  q->m1=c*(1.-c-d);
  q->n0=-c*d;
  q->n1=d*(1.-d);
  q->n2=c*(1.-c);

  q->t_base  =  tref;
  q->p_base  =  pref;

//   q->DP[0] =  (double) R *  (t2-t3)*d2r;
//   q->DP[1] =  (double) R *  (t3-t1)*d2r;
//   q->DP[2] =  (double) R *  (t1-t2)*d2r;
//   q->DQ[0] =  (double) R *  (p2-p3)*d2r;
//   q->DQ[1] =  (double) R *  (p3-p1)*d2r;
//   q->DQ[2] =  (double) R *  (p1-p2)*d2r;
/*------------------------------------------------------------------------------
  WARNING, PSEUDO area (actual area = pseudo*cos(latitude))
  note: PSEUDO Area=0.5*ef_jacobien*R²*dtr² */
// The area of a quadrilateral ABCD can be calculated using vectors. Let vectors AC and BD form the diagonals from A to C and from B to D. The area of the quadrilateral is then
//
//     K = \tfrac{1}{2} |\mathbf{AC}\times\mathbf{BD}|,
  u=vector2D_t(t1,p1,t4,p4);
  v=vector2D_t(t2,p2,t3,p3);
  q->Area  =  0.5*abs(u*v)*d2r*d2r*R*R;

//  q->Area = 0.5*R* (t1*q->DQ[0] + t2*q->DQ[1] + t3*q->DQ[2])*d2r;

  if(isnan(q->Area)==1) {
    q->Area=0.;
    }
  
  q->cosinus=cos((p1+p2+p3+p4)*d2r/4.);

  q->TrueArea =q->Area*q->cosinus;
  
//   a1=t2-t1;
//   a2=t1-t2+t3-t4;
//
//   b1=t4-t1;
//   b2=t1-t2+t3-t4;
//
//   c1=p2-p1;
//   c2=p1-p2+p3-p4;
//
//   d1=p4-p1;
//   d2=p1-p2+p3-p4;
  
  a1=t2-t1;
  a2=t1-t2+t4-t3;
  
  b1=t3-t1;
  b2=t1-t2+t4-t3;
  
  c1=p2-p1;
  c2=p1-p2+p4-p3;
  
  d1=p3-p1;
  d2=p1-p2+p4-p3;
  
  double test=(a1*d1 - b1*c1) + 0.5*(a1*d2 -c1*b2 + a2*d1 - b1*c2) + 0.25*(a2*d2 - b2*c2);
  test *=d2r*d2r*R*R;
  
  
  q->J[0]=a1*d1 - b1*c1;
  
  q->J[1]=(a1*d1 - b1*c1) + (a1*d2 -c1*b2);
  
  q->J[2]=(a1*d1 - b1*c1) + (a1*d2 -c1*b2) + (a2*d1 - b1*c2) + (a2*d2 - b2*c2);
  
  q->J[3]=(a1*d1 - b1*c1) + (a2*d1 - b1*c2);
  
  for(int k=0;k<4;k++) q->J[k]*=d2r*d2r*R*R;
  
  double p[4]={1,1,1,1};
  test=fe_integrale_QP1xQP1_2D(p,q->J);
//  test *=d2r*d2r*R*R;
  
//   for(i=0;i<3;i++) {
//     n1=q->vertex[i];
//     n2=q->vertex[(i+1)%3];
//     t1=mesh.vertices[n1].lon;
//     t2=mesh.vertices[n2].lon;
//
//     if(mesh.type==0) t2=geo_recale(t2,t1,(double) 180.0);
//
//     p1=mesh.vertices[n1].lat;
//     p2=mesh.vertices[n2].lat;
//     alpha = (double) 0.5 * (p1+p2)*d2r;
//     dlon=t2-t1;
//     dlat=p2-p1;
//     dX = R * cos(alpha) * dlon*d2r;
//     dY = R * dlat*d2r;
//     ds = sqrt(dX * dX + dY * dY);
// /*----------------------------------------------------------------------
//     outward normal is clockwise rotated */
//     q->nx[i] = +dY / ds;
//     q->ny[i] = -dX / ds;
//     q->l[i] = ds;
//     }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_test_quadrangle_area(const mesh_t & mesh,int m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(mesh.quadrangles[m].Area==0.0) {
    TRAP_ERR_RETURN(MESH_STATUS_FLAT_ELEMENT,1,"%s:quadrangle %d is flat !\n",__func__,m);
    }
  if(isnan(mesh.quadrangles[m].Area)) {
    TRAP_ERR_RETURN(MESH_STATUS_FLAT_ELEMENT,1,"%s:quadrangle %d is WRONG !\n",__func__,m);
    }
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_intpl_QLP1_template(const mesh_t & mesh,const T *buffer, T mask, double t, double p, int m, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T zz;
  double x, y, beta[4];
  int nodes[4];
  int i;
  int status;

  status=fe_test_quadrangle_area(mesh,m);
  if(status)TRAP_ERR_RETURN(status,1,"%s:fe_test_quadrangle_area(,%d)=%d\n",__func__,m,status);
  status = fe_affine_directe(mesh.quadrangles[m],t, p, &x, &y, mesh.type);

  fe_QP1base(x, y, beta);

  if(status == 0) {
    for(i = 0; i < 4; i++) {
      nodes[i] = mesh.quadrangles[m].vertex[i];
      }
    }
  else {
    printf("error in ef_beta: %f %f %f %f\n", t, p, x, y);
    }

  zz = 0;
  for(i = 0; i < 4; i++) {
    zz += buffer[nodes[i]] * (T) beta[i];
    }

  *z = zz;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_intpl_QLP1(const mesh_t & mesh,const float *buffer,float mask,double t,double p,int m,float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_QLP1_template(mesh, buffer, mask, t, p, m, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_intpl_QLP1(const mesh_t & mesh,const double *buffer,double mask,double t,double p,int m,double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_QLP1_template(mesh, buffer, mask, t, p, m, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_intpl_QLP1(const mesh_t & mesh,const complex<float> *buffer,complex<float> mask,double t,double p,int m,complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_QLP1_template(mesh, buffer, mask, t, p, m, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_intpl_QLP1(const mesh_t & mesh,const complex<double> *buffer,complex<double> mask,double t,double p,int m,complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_QLP1_template(mesh, buffer, mask, t, p, m, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_intpl_CQN1_template(const mesh_t & mesh,const T *buffer, T mask, double t, double p, int m, T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T zz;
  double x, y, beta[4];
  int nodes[4];
  int i;
  int status;

  if(mesh.quadrangles[m].Area==0.0) {
    status = fe_initaffine_spherical(mesh, &(mesh.quadrangles[m]),m);
    if(status != 0) goto error;
    status = fe_affine_directe(mesh.quadrangles[m],t, p, &x, &y, mesh.type);
    }
  else {
    status = fe_affine_directe(mesh.quadrangles[m],t, p, &x, &y, mesh.type);
    }

  fe_QP1base(x, y, beta);

  if(status == 0) {
    for(i = 0; i < 4; i++) {
      nodes[i] = mesh.quadrangles[m].edges[i];
      }
    }
  else {
    printf("error in ef_beta: %f %f %f %f\n", t, p, x, y);
    }

  zz = 0;
  for(i = 0; i < 4; i++) {
    int n1=nodes[i];
    int n2=nodes[(i+3)%4];
    zz += (T) 0.5*(buffer[n1]+buffer[n2]) * (T) beta[i];
    }

  *z = zz;
  return (status);

error:
  printf("error in fe_intpl_CQN1() \n");
  return (-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_intpl_CQN1(const mesh_t & mesh,const float *buffer, float mask, double t, double p, int m, float *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_CQN1_template(mesh, buffer, mask, t, p, m, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_intpl_CQN1(const mesh_t & mesh,const double *buffer, double mask, double t, double p, int m, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_CQN1_template(mesh, buffer, mask, t, p, m, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_intpl_CQN1(const mesh_t & mesh,const complex<float> *buffer, complex<float> mask, double t, double p, int m, complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_CQN1_template(mesh, buffer, mask, t, p, m, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_intpl_CQN1(const mesh_t & mesh,const complex<double> *buffer, complex<double> mask, double t, double p, int m, complex<double> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=fe_intpl_CQN1_template(mesh, buffer, mask, t, p, m, z);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_e2n (mesh_t *mesh,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  build the vertices' neighbours list from the elements list

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i,j,k,m,n,n1,maxn;
  
  struct timeval b4;
  gettimeofday(&b4);
  
/* NOTE (Damien Allain, 2013-08-26)
  Doing things thrice is actually faster than when using e.g. vector<int>'s... */

/*------------------------------------------------------------------------------
  rough estimate of the maximum count of neighbours*/
  for(n=0;n<mesh->nvtxs;n++) mesh->vertices[n].nngh=0;
  if(verbose>1) STDERR_BASE_LINE("%gs\n",difftime(&b4));

  for(m=0;m<mesh->ntriangles;m++)
    for(i=0;i<3;i++) {
      n=mesh->triangles[m].vertex[i];
      mesh->vertices[n].nngh+=2;
      }
  if(verbose>1) STDERR_BASE_LINE("%gs\n",difftime(&b4));

  for(m=0;m<mesh->nquadrangles;m++)
    for(i=0;i<4;i++) {
      n=mesh->quadrangles[m].vertex[i];
      mesh->vertices[n].nngh+=2;
      }
  if(verbose>1) STDERR_BASE_LINE("%gs\n",difftime(&b4));

/*------------------------------------------------------------------------------
  precise estimate of the maximum count of neighbours*/
  for (n=0; n<mesh->nvtxs; n++) {
    vertex_t *vertexn=&mesh->vertices[n];
    vertexn->ngh=aset(vertexn->nngh,-1);
    vertexn->nngh=0;
    }
  if(verbose>1) STDERR_BASE_LINE("%gs\n",difftime(&b4));
  
  for(m=0;m<mesh->ntriangles;m++) {
    const triangle_t *trianglem=&mesh->triangles[m];
    for(i=0;i<3;i++) {
      n=trianglem->vertex[i];
      vertex_t *vertexn=&mesh->vertices[n];
      for(j=0;j<3;j++) {
        if(i==j) continue;
        n1=trianglem->vertex[j];
        k=pos(n1,vertexn->ngh,vertexn->nngh);
        if(k==-1) {
          vertexn->ngh[vertexn->nngh]=n1;
          vertexn->nngh++;
          }
        }
      }
    }
  if(verbose>1) STDERR_BASE_LINE("%gs\n",difftime(&b4));

  for(m=0;m<mesh->nquadrangles;m++) {
    for(i=0;i<4;i++) {
      n=mesh->quadrangles[m].vertex[i];
      for(j=0;j<4;j++) {
        if(abs(i-j)%2==0) continue;
        n1=mesh->quadrangles[m].vertex[j];
        k=pos(n1,mesh->vertices[n].ngh,mesh->vertices[n].nngh);
        if(k==-1) {
          mesh->vertices[n].ngh[mesh->vertices[n].nngh]=n1;
          mesh->vertices[n].nngh++;
          }
        }
      }
    }
  if(verbose>1) STDERR_BASE_LINE("%gs\n",difftime(&b4));

/*------------------------------------------------------------------------------
  shorten allocation of vertices' neighbours list */
  maxn=0;
  for (n=0; n<mesh->nvtxs; n++) {
    vertex_t *vertexn=&mesh->vertices[n];
    
    int *oldNgh=vertexn->ngh;
    vertexn->ngh=new int[vertexn->nngh];
    valcpy(vertexn->ngh,oldNgh,vertexn->nngh);
    delete[] oldNgh;
    
    updatemax(&maxn, vertexn->nngh);
    }
  if(verbose>1) STDERR_BASE_LINE("%gs\n",difftime(&b4));
  
  mesh->nnghm=maxn;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_edgetable_Q(mesh_t *mesh, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   k,l,j,j1,j2,n1,n2,n,status;
  int   nedges,chk=0,count=0;
  int   nelt,nndes;
  int   interior=0,boundary=0,weird=0,orphans=0;
  int   *ncells,ncellsmax=0,**cells;
  double t1,t2,p1,p2;
  edge_t *edges;

retry:
  nedges=0;

  nelt=mesh->nelements;
  nndes=mesh->nvtxs;

/*-----------------------------------------------------------------------------
  count edges*/
  for(n=0; n<nndes; n++) {
    for(j=0; j<mesh->vertices[n].nngh;j++) {
      n2=mesh->vertices[n].ngh[j];
      if(n2>n) nedges++;
      }
    }

/*-----------------------------------------------------------------------------
  count max number of element connected to a given vertex*/
  ncells=new int[nndes];
  cells=new int*[nndes];

  for(n=0; n<nndes; n++) ncells[n]=0;

  for(k=0; k<mesh->ntriangles; k++) {
    for(j=0; j<3;j++) {
      n=mesh->triangles[k].vertex[j];
      ncells[n]++;
      }
    }
  for(k=0; k<mesh->nquadrangles; k++) {
    for(j=0; j<4;j++) {
      n=mesh->quadrangles[k].vertex[j];
      ncells[n]++;
      }
    }

  for(n=0; n<nndes; n++) {
    if(ncells[n] ==0) {
      printf("%s: unused node %d %d\n",__func__,n,mesh->vertices[n].nngh);
      }
    }

  for(n=0; n<nndes; n++) updatemax(&ncellsmax,ncells[n]);
  for(n=0; n<nndes; n++) cells[n]=new int[ncellsmax];
  for(n=0; n<nndes; n++) ncells[n]=0;

/*-----------------------------------------------------------------------------
  build element connection table*/
  for(k=0; k<mesh->ntriangles; k++) {
    for(j=0; j<3;j++) {
      n=mesh->triangles[k].vertex[j];
      cells[n][ncells[n]]=k;
      ncells[n]++;
      }
    }
  for(k=0; k<mesh->nquadrangles; k++) {
    for(j=0; j<4;j++) {
      n=mesh->quadrangles[k].vertex[j];
/* *-----------------------------------------------------------------------------
      implicitely (and dangerously) take element list in triangles then quadrangles order*/
      cells[n][ncells[n]]=k+mesh->ntriangles;
      ncells[n]++;
      }
    }

  edges=new edge_t[nedges];
  count=0;
  for (n=0;n<nedges;n++) {
    edges[n].nshared=0;
    edges[n].shared[0]=-1;
    edges[n].shared[1]=-1;
    }

  for(n1=0; n1<nndes; n1++) {
    mesh->vertices[n1].nelmts=ncells[n1];
    mesh->vertices[n1].elmts=cells[n1];
    for(j=0; j<mesh->vertices[n1].nngh;j++) {
      n2=mesh->vertices[n1].ngh[j];
      if(n2>n1) {
/*-----------------------------------------------------------------------------
        create a new edge entry*/
        edges[count].extremity[0]=n1;
        edges[count].extremity[1]=n2;
/*-----------------------------------------------------------------------------
        scan elements shared by n1 and n2*/
        for (k=0; k<ncells[n1];k++)
          for (l=0; l<ncells[n2];l++)
            if(cells[n1][k]==cells[n2][l]) {
              if (edges[count].nshared >1) {
                printf ("found edge %d to be weird: %d %d, element %d\n",count,n1,n2,cells[n1][k]);
                printf ("already in list: elements %d %d\n",edges[count].shared[0],edges[count].shared[1]);
                }
              else {
                edges[count].shared[edges[count].nshared]=cells[n1][k];
                edges[count].nshared++;
                }
              }
        count++;
        }
      }
    }
  
  /* NOT USING deletep2D(&cells,) BECAUSE ABOVE DOES mesh->vertices[n1].elmts=cells[n1] */
  deletep(&cells);
  deletep(&ncells);

  boundary=0;
  interior=0;
  weird=0;
  orphans=0;
  for (n=0;n<nedges;n++) {
    if (edges[n].nshared==0) {
      weird++;
      n1=edges[n].extremity[0];
      n2=edges[n].extremity[1];
      if(verbose) {
        printf ("found %d weird edges, no element attached: edge=%d vertex=%d vertex=%d\n",weird,n,n1,n2);
        printf ("position: %lf %lf\n",mesh->vertices[n1].lon,mesh->vertices[n1].lat);
        }
      status=fe_disconnectvertices(*mesh, n1, n2);
      orphans++;
      }
    if (edges[n].nshared >2) {
      weird++;
      n1=edges[n].extremity[0];
      n2=edges[n].extremity[1];
      if(verbose) {
        printf ("found %d weird edges, too many elements attached: edge=%d vertex=%d vertex=%d\n",weird,n,n1,n2);
        printf ("position: %lf %lf\n",mesh->vertices[n1].lon,mesh->vertices[n1].lat);
        }
      }
    if (edges[n].nshared==1) boundary++;
    if (edges[n].nshared==2) interior++;
    }
  
  if(orphans!=0) {
    edges->destroy();
    if(verbose) printf ("found %d orphans, vertices neighbours cleaned, re-build edge set\n", orphans);
    goto retry;
    }
    
/**-----------------------------------------------------------------------------
  
  relationship between elements and edges:
  
  2*interior+boundary = nquadrangles*4 + ntriangles*3
                      = nquadrangles   + nelts*3

-----------------------------------------------------------------------------*/

  if((verbose==1) || (weird!=0)) {
    STDOUT_BASE_LINE_FUNC("\nfound %d total edges\n",nedges);
    printf ("found %d interior edges\n",interior);
    printf ("found %d boundary edges\n",boundary);
    printf ("found %d weird edges (should be zero) \n",weird);
    printf ("found presumedly %d edges, predicted %d\n",(interior*2+boundary),(mesh->ntriangles*3+mesh->nquadrangles*4));
    }

//  if(weird!=0) return(-1);
  
/*-----------------------------------------------------------------------------
  scan boundary edges for outward normal determination*/
  for (n=0;n<nedges;n++) {
    n1=edges[n].extremity[0];
    n2=edges[n].extremity[1];
/*-----------------------------------------------------------------------------
    create a new node at the middle of the edge */
    t1=mesh->vertices[n1].lon;
    t2=mesh->vertices[n2].lon;
    if(mesh->type==0) t2=geo_recale(t2,t1,(double) 180.0);
    p1=mesh->vertices[n1].lat;
    p2=mesh->vertices[n2].lat;
    edges[n].lon=0.5*(t1+t2);
    edges[n].lat=0.5*(p1+p2);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  scan boundary edges for outward normal determination
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(n = 0; n < nedges; n++) {
    if(edges[n].nshared == 1) {
      n1 = edges[n].extremity[0];
      n2 = edges[n].extremity[1];
      for(j1 = 0; j1 < 4; j1++)
        if(mesh->quadrangles[edges[n].shared[0]].vertex[j1] == n1)
          break;
      for(j2 = 0; j2 < 4; j2++)
        if(mesh->quadrangles[edges[n].shared[0]].vertex[j2] == n2)
          break;
      if((j2 - j1 == -1) or (j2 - j1 == 3)) {
/*-----------------------------------------------------------------------------
        put in direct rotation, so normal vector will be easier to compute*/
        edges[n].extremity[0] = n2;
        edges[n].extremity[1] = n1;
        }
      }
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  associate appropriate edges to each element (kind of reverse table) 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(k=0; k<nelt; k++) {
    for(j=0; j<4;j++) {
      mesh->quadrangles[k].edges[j]=-1;
      }
    }

  for (n=0;n<nedges;n++) {
    for (k=0;k<edges[n].nshared;k++) {
/*       ptr=edges[n].shared[k]; */
//      ptr=&(elt[edges[n].shared[k]]);

      int m =edges[n].shared[k];
      n1=edges[n].extremity[0];
      n2=edges[n].extremity[1];
      
      j1=vpos(n1, mesh->quadrangles[m].vertex, 4);
      j2=vpos(n2, mesh->quadrangles[m].vertex, 4);

/**-----------------------------------------------------------------------------
      node 0-1 : side 0, but edge 2 (opposite to node 2)
      node 1-2 : side 1, but edge 0
      node 2-0 : side 2, but edge 1*/
      if((j2-j1==1)||(j2-j1==-3)) {
/*-----------------------------------------------------------------------------
        direct orientation, use j1 (smaller node index in element) to switch edge index*/
        switch(j1) {
          case 0:
            mesh->quadrangles[m].edges[0]=n;
            break;
          case 1:
            mesh->quadrangles[m].edges[1]=n;
            break;
         case 2:
            mesh->quadrangles[m].edges[2]=n;
            break;
         case 3:
            mesh->quadrangles[m].edges[3]=n;
            break;
          }
        edges[n].Tx= -mesh->quadrangles[m].ny[j1]; /*warning initialisation not guaranted*/
        edges[n].Ty=  mesh->quadrangles[m].nx[j1];
        edges[n].L =  mesh->quadrangles[m].l [j1];
        }
      else {
/*-----------------------------------------------------------------------------
        reverse orientation, use j2 (smaller node index in element) to switch edge index*/
        switch(j2) {
          case 0:
            mesh->quadrangles[m].edges[0]=n;
            break;
          case 1:
            mesh->quadrangles[m].edges[1]=n;
            break;
         case 2:
            mesh->quadrangles[m].edges[2]=n;
            break;
         case 3:
            mesh->quadrangles[m].edges[3]=n;
            break;
          }
        edges[n].Tx=  mesh->quadrangles[m].ny[j2];
        edges[n].Ty= -mesh->quadrangles[m].nx[j2];
        edges[n].L =  mesh->quadrangles[m].l [j2];
        }
      }
    }

/*-----------------------------------------------------------------------------
  cheks if previous step is correct */
  chk=0;
  for(k=0; k<nelt; k++) {
    for(j=0; j<4;j++) {
      if(mesh->quadrangles[k].edges[j]==-1) {
        chk++;
        printf ("found anomaly %d at element %d, edge %d\n",chk,k,j);
        }
      }
    }

//   delete[] ncells;
//   for(n=0; n<nndes; n++) delete[] cells[n];
//   delete[] cells;
 
  mesh->nedges=nedges;
  mesh->edges=edges;

//   status = fe_vertex_crosstables02(mesh); /// HERE !!!
  status=fe_Vertex_EdgesXtables(mesh);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_SeekQuadrangles(mesh_t *mesh, Polygons::ElementaryCycles &elementaryCycles,int maxsize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Make an element list from nodes' neighbours description.
  If consecutive neighbours of a node are connected, they
  form an element with the node.
----------------------------------------------------------------------*/
{
  int i,j,k,l,m,n;
  int rotation=-1;
  int max_nb,nndes,*flag,*type,*size;
  int ntriangles,nquadrangles,npentacles,nelements;
  int **tmp;
  int status;
  int nnpe[100];
  int nodes[maxsize];
  int cardinal[maxsize];

  if(mesh->type==-1) mesh->type=0;

  nndes=mesh->nvtxs;
  max_nb=mesh->nnghm;

//   int nprocs=initialize_OPENMP(-1);
// #ifdef OPEN_MP
//   #pragma omp parallel for private(status) if(nprocs>1)
// #endif
//   for(i=0; i<nndes; i++) {
//     status=fe_order(*mesh,i,(double) -1., (double) 0.,rotation);
//     }

  flag=new int[nndes];
  for(i=0; i<nndes; i++) flag[i]=0;

  for(i=0; i<nndes; i++) {
    if(mesh->vertices[i].nngh < 2 ) {
      printf("fe_SeekQuadrangles: pathetic node %d\n",i);
      flag[i]=1;
      }
    }

  tmp= new int*[2*nndes];
  type=new int[2*nndes];
  size=new int[2*nndes];
  
/**-----------------------------------------------------------------------------
  this makes the algorithm unconvenient for polygone detection */  
  for(i=0; i<2*nndes; i++) tmp[i]=new int[maxsize-1];

  ntriangles   = 0;
  nquadrangles = 0;
  npentacles   = 0;
  nelements    = 0;

  for(l=0;l<maxsize;l++) cardinal[l]=0;

/**-----------------------------------------------------------------------------
  seek anti-clockwise cycles with limited maxise */  
  for(n=0;n<nndes;n++) {
    if(flag[n]!=0) continue;
    status=fe_order(*mesh,n,(double) -1., (double) 0.,rotation);  // not necessary
    nodes[0]=n;
    for(k=0; k<mesh->vertices[n].nngh; k++) {
      nodes[1]=mesh->vertices[nodes[0]].ngh[k];
/**-----------------------------------------------------------------------------
      to eliminate redundant cycle, keep the one starting with the lowest node index */  
      if(nodes[1] <= nodes[0]) continue;
      status = fe_order2(*mesh, nodes[1], nodes[0], rotation);
      nodes[2]=mesh->vertices[nodes[1]].ngh[0];
      if(nodes[2] == nodes[1]) continue;
      if(nodes[2] <= nodes[0]) continue;
      for(l=3;l<maxsize;l++) {
        status = fe_order2(*mesh, nodes[l-1], nodes[l-2], rotation);
        nodes[l]=mesh->vertices[nodes[l-1]].ngh[0];
        if(nodes[l] < nodes[0]) goto skip;
/* *----------------------------------------------------------------------------
        to avoid mis-formed cycles*/
        if(nodes[l]==nodes[0]) {
          double angle, sum=0;
          for(j=0;j<l;j++) {
            angle=fe_angle(*mesh, nodes[(j+1) % l], nodes[j], nodes[(j+2) % l]);
            sum+=angle;
            }
          if(sum>0) {
            goto skip;
            }
          for(j=0;j<l;j++) tmp[nelements][j]=nodes[j];
          size[nelements]=l;
          cardinal[l]++;
          nelements++;
          goto skip;
          }
        }
skip:
      continue;
      }
    }

  for(m=0;m<nelements;m++) {
    switch (size[m]){
      case 3:
        type[m]=FE_TRIANGLE;
        ntriangles++;
        break;
      case 4:
        type[m]=FE_QUADRANGLE;
        nquadrangles++;
        break;
      case 5:
        type[m]=FE_PENTACLE;
        npentacles++;
        break;
      default:
        type[m]=FE_UNDEFINED;
        break;
      }
    }
  nnpe[FE_TRIANGLE]  =3;
  nnpe[FE_QUADRANGLE]=4;
  nnpe[FE_PENTACLE]  =5;

  for(i=0; i<nndes; i++) flag[i]=-1;

  for(m=0;m<nelements;m++) {
    Polygons::Adjacency cycle;
    for(k=0;k<size[m];k++) {
      flag[tmp[m][k]]=0;
      cycle.push_back(tmp[m][k]);
      }
    elementaryCycles.push_back(cycle);
    }

  for(i=0; i<nndes; i++) {
    if(flag[i]==-1) {
      STDOUT_BASE_LINE_FUNC("unused node %d\n",i);
      }
    }

  for(m=0;m<nelements;m++) {
    if(size[m]>5) {
      printf("oversized element %d: 1st node %d\n",m,tmp[m][0]);
      }
    }
    
  for(i=0; i<2*nndes; i++) delete[] tmp[i];
  delete[] tmp;
  delete[] flag;
  delete[] size;
  delete[] type;

  return(0);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int quadrangle_readneigh (const string &polygonsFileName, mesh_t & final) throw ()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
using namespace Polygons;

  plg_array_t set;
  int k,m,n,status;
  size_t add_ntriangles,add_nquadrangles;
  ElementaryCycles elementaryCycles;
//  gmesh_t<quadrangle_t> mesh;
  mesh_t mesh;
  mesh_t tmp;
  size_t size,count;
//  Adjacency cycle;

  try {
/* *----------------------------------------------------------------------------
    Read graph from polygonsFileName */

    status=fe_readmesh(polygonsFileName.c_str(), MESH_FILE_FORMAT_TRIGRID, &tmp);

    status=fe_SeekQuadrangles(&tmp, elementaryCycles, 6);

    GraphBase<double> graph(polygonsFileName);
    printf("found %d nodes \n",graph.size());

/* *----------------------------------------------------------------------------
    search cycles from graph*/
//    PolygonSet polygons = searchPolygons<double,double>(graph,5);
//    elementaryCycles = searchCycles<double>(graph,5);

    mesh.nquadrangles=0;
    mesh.ntriangles=0;
    mesh.nelements=0;
    add_ntriangles=0;
    add_nquadrangles=0;
    for (ElementaryCycles::const_iterator i = elementaryCycles.begin(); i != elementaryCycles.end(); i++) {
      // The current cycle.
      const Adjacency cycle = *i;
//      cycle = *i;
      size=cycle.size();
      switch(size) {
        case 0:
        case 1:
        case 2:
          printf("Undersized cycle %d\n",size);
          break;
        case 3:
          mesh.ntriangles++;
          mesh.nelements++;
          break;
        case 4:
          mesh.nquadrangles++;
          mesh.nelements++;
          break;
        case 5:
          add_ntriangles++;
          mesh.nelements++;
          add_nquadrangles++;
          mesh.nelements++;
          break;
        default:
          printf("Oversized cycle %d\n",size);
          break;
        }
      }

/* *----------------------------------------------------------------------------
    create vertices (geometrical nodes)*/
//     mesh.nvtxs=graph.size();
//     mesh.vertices=new vertex_t[mesh.nvtxs];
//     for(n=0;n<mesh.nvtxs;n++) {
//       PointBase<double> p;
//       p=graph[n].coordinates();
//       mesh.vertices[n].lon=p.x();
//       mesh.vertices[n].lat=p.y();
//       }

    mesh.nvtxs=tmp.nvtxs;
    mesh.vertices=new vertex_t[mesh.nvtxs];
    for(n=0;n<mesh.nvtxs;n++) {
      mesh.vertices[n].lon=tmp.vertices[n].lon;
      mesh.vertices[n].lat=tmp.vertices[n].lat;
      }

/* *----------------------------------------------------------------------------
    create triangles*/
    mesh.triangles=new triangle_t[mesh.ntriangles+add_ntriangles];

    count=0;
    for (ElementaryCycles::const_iterator i = elementaryCycles.begin(); i != elementaryCycles.end(); i++) {
      // The current cycle.
      const Adjacency cycle = *i;
      // place is a reference to a node of the cycle.
      size=cycle.size();
      if(size==3) {
        k=0;
        for (Adjacency::const_iterator place = cycle.begin(); place != cycle.end(); place++) {
          mesh.triangles[count].vertex[k]=*place;
          k++;
          }
        count++;
        }
      }

/* *----------------------------------------------------------------------------
    create quadrangles*/
    mesh.quadrangles=new quadrangle_t[mesh.nquadrangles+add_nquadrangles];

    count=0;
    for (ElementaryCycles::const_iterator i = elementaryCycles.begin(); i != elementaryCycles.end(); i++) {
      // The current cycle.
      const Adjacency cycle = *i;
      // place is a reference to a node of the cycle.
      size=cycle.size();
      if(size==4) {
        k=0;
        for (Adjacency::const_iterator place = cycle.begin(); place != cycle.end(); place++) {
          mesh.quadrangles[count].vertex[k]=*place;
          k++;
          }
        count++;
        }
      }

/* *----------------------------------------------------------------------------
    split pentacles into 1 quadrangle and 1 triangle */
    count=0;
    for (ElementaryCycles::const_iterator i = elementaryCycles.begin(); i != elementaryCycles.end(); i++) {
      // The current cycle.
      const Adjacency cycle = *i;
      int j,found=0,start;
      // place is a reference to a node of the cycle.
      size=cycle.size();
      if(size==5) {
        Adjacency::const_iterator place,last=cycle.end();
        last--;
//        Adjacency tmp;
/* *----------------------------------------------------------------------------
        find pentacle's head */
        k=0;
        for (place = cycle.begin(); place != cycle.end(); place++) {
          n=*place;
          Adjacency p=graph[n].adjacency();
          if(p.size()==4) {
            found=1;
            start=k;
            break;
            }
          k++;
          }
        if(found==0) {
          printf("pentacle (cycle %d) not a connexion element\n",count);
          for (Adjacency::const_iterator place = cycle.begin(); place != cycle.end(); place++) {
            printf("%d ",*place);
            }
          printf("\n");
          start=0;
          }
/* *----------------------------------------------------------------------------
        i.e. start-1 */
//           for (Adjacency::const_iterator place = cycle.begin(); place != cycle.end(); place++) {
//             printf("%d ",*place);
//             }
//           printf("\n");
        k=(start+4)%5;
        place=cycle.begin();
        for (j=0; j<k; j++) {place++;};
        for (j=0; j<3; j++) {
          mesh.triangles[mesh.ntriangles].vertex[j]=*place;
          if(place==last) {
            place=cycle.begin();
            }
          else {
            place++;
            }
          }
        mesh.ntriangles++;
/* *----------------------------------------------------------------------------
        i.e. start+1 */
        k=(start+1)%5;
        place=cycle.begin();
        for (j=0; j<k; j++) {place++;};
        for (j=0; j<4; j++) {
          mesh.quadrangles[mesh.nquadrangles].vertex[j]=*place;
          if(place==last) place=cycle.begin();
          else place++;
          }
        mesh.nquadrangles++;
        count++;
        }
      }
    }

/* *----------------------------------------------------------------------------
  Part of the code executed if a problem occur during file reading*/
  catch (TugoExceptions::ReadError file) {
    TRAP_ERR_EXIT(-2, "An error occured when reading the file \"%s\".",file.fileName().c_str());
    }

/* *----------------------------------------------------------------------------
  testing interpolation*/
  double *buffer=new double[mesh.nvtxs],mask=1.e+10;
  for(n=0;n<mesh.nvtxs;n++) {
    buffer[n]=mesh.vertices[n].lon;
    }
  for(m=0;m<mesh.nquadrangles;m++) {
    status=fe_initaffine_spherical(mesh, &(mesh.quadrangles[m]), m);
    double t=3.7929, p= 43.2974,check;
    double ksi,eta;
    t=mesh.vertices[3].lon;
    p=mesh.vertices[3].lat;
    int detect = fe_testelement_v2(mesh,mesh.quadrangles[m], t, p,0);
    if(detect==1) {
      status= fe_affine_directe(mesh.quadrangles[m], t, p, &ksi, &eta,0);
      status=fe_intpl_QLP1( mesh, buffer, mask, t, p, m, &check);
      printf("detection/interpolation test: quadrangle %d detected, ksi=%lf eta=%lf check=%lf\n",m,ksi,eta, check);
      }
    }
  delete[] buffer;

  status=fe_e2n (&mesh);
  status=fe_savemesh("test.nei",MESH_FILE_FORMAT_TRIGRID,tmp);
  status=fe_savemesh("processed.nei",MESH_FILE_FORMAT_TRIGRID,mesh);
  status=fe_edgetable_Q(&mesh,0);

  mesh.nlayers=1;
  mesh.nlevels=2;
  status=fe_savemeshNC2D_new("test.nc",  mesh, 1);
  
  final=mesh;

// 3.79 43.32
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_readQmesh_ascii(const string &polygonsFileName, mesh_t *mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int count;
  plg_array_t set;
  mesh_t tmp,final;
  string cleaned;
  size_t pos;
  const string& str(".cleaned");

/* *----------------------------------------------------------------------------
  get polygons from neighbour file*/
//  set= plg_readneigh (polygonsFileName);
  status=fe_readmesh(polygonsFileName.c_str(),MESH_FILE_FORMAT_TRIGRID,&tmp);

/* *----------------------------------------------------------------------------
  force neigbours list's consistency (neighbourship symetry)*/
  count=fe_chkvertex_consistency(&tmp,1);
  count=fe_chkvertex_consistency(&tmp,1);
 
/* *----------------------------------------------------------------------------
  remove vertices with no connexions at all*/
  status=fe_cleanvertices(&tmp, false);

/* *----------------------------------------------------------------------------
  save "consistent" mesh*/
  cleaned=polygonsFileName;
  pos=cleaned.rfind(string(".nei"));
  cleaned.insert(pos,str);
  status=fe_savemesh(cleaned.c_str(),MESH_FILE_FORMAT_TRIGRID,tmp);

/* *----------------------------------------------------------------------------
  now process cycles*/
  status=quadrangle_readneigh(cleaned.c_str(),final);

/* *----------------------------------------------------------------------------
  convert polygons into mesh*/
//  mesh->nelts=set.n;
  return(status);
}

