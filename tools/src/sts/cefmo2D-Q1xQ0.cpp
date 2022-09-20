
/**************************************************************************

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

***************************************************************************/

#include "tugo-prototypes.h"
#include "cefmo.h"
#include "matrix.h"
#include "tides.h"

#define Z_NNPE 1


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpGradient2D_CQN1xCQP0_02(const mesh_t & mesh, cefmo_t & cefmo, parameter_t & data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 edge x mid-point version

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

{
  int i, j, k, l, m, mm, n, row;
  int status;
  int udim,zdim;
  double S, Lmm, Cd, h;
  double p[Z_NNPE],sumx,sumy,C;
  ordering_t ordering;
  int *pointer,*incidence,*cardinal;
  double dzdx,dzdy,factor;
  quadrangle_t quadrangle;

  discretisation_t & u_descriptor=*cefmo.u_descriptor;
  discretisation_t & z_descriptor=*cefmo.z_descriptor;

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);

  pointer   = cefmo.G1.ordering->pointer;
  incidence = cefmo.G1.ordering->incidence;
  cardinal  = cefmo.G1.ordering->cardinal;

  cefmo.u_LumpedMW=new double[udim];
  for(n=0;n<udim;n++) cefmo.u_LumpedMW[n]=0.0;
  
/**-----------------------------------------------------------------------
  compute pressure gradients matrices G1 & G2*/
  for(m = 0; m < mesh.nquadrangles; m++) {
    quadrangle=mesh.quadrangles[m];
/**----------------------------------------------------------------------
    transport version */
    h=cefmo.state.h_z[m];
    for(i = 0; i < u_descriptor.nnpe; i++) {
      n = mesh.quadrangles[m].edges[i];
      cefmo.u_LumpedMW[n]+=0.5*quadrangle.Area;
      if(mesh.edges[n].nshared==1) continue;
/**----------------------------------------------------------------------
      element summation, half-discontinuity */
      sumx =  0.5*h*quadrangle.l[i]*quadrangle.nx[i];
      sumy =  0.5*h*quadrangle.l[i]*quadrangle.ny[i];
      
      size_t start=pointer[n];
      row=vpos(m, &(incidence[start]), cardinal[n]);
      if(row==-1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
      cefmo.G1.packed[pointer[n]+row]-= P_g *sumx;
      cefmo.G2.packed[pointer[n]+row]-= P_g *sumy;
      if(mesh.edges[n].shared[0]==m) mm=mesh.edges[n].shared[1];
      else mm=mesh.edges[n].shared[0];
      row=vpos(mm, &(incidence[start]), cardinal[n]);
      if(row==-1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
      cefmo.G1.packed[pointer[n]+row]+= P_g *sumx;
      cefmo.G2.packed[pointer[n]+row]+= P_g *sumy;

/**----------------------------------------------------------------------------
      half-side-length (0.5), but double weight (1.0) for open boundary neighbours*/
      n=  mesh.quadrangles[m].edges[(i+2)%4]; /// opposite edge
      if(mesh.edges[n].nshared==1) factor=1;
      else factor=0.5;
      
      n = mesh.quadrangles[m].edges[(i+1)%4];
//      if(mesh.edges[n].nshared!=1) {
      start=pointer[n];
      row=vpos(m, &(incidence[start]), cardinal[n]);
      if(row==-1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
      cefmo.G1.packed[pointer[n]+row]-= factor*P_g *sumx;
      cefmo.G2.packed[pointer[n]+row]-= factor*P_g *sumy;
      row=vpos(mm, &(incidence[start]), cardinal[n]);
      if(row==-1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
      cefmo.G1.packed[pointer[n]+row]+= factor*P_g *sumx;
      cefmo.G2.packed[pointer[n]+row]+= factor*P_g *sumy;
//      }
      
      n = mesh.quadrangles[m].edges[(i+3)%4];
//      if(mesh.edges[n].nshared!=1) {
      start=pointer[n];
      row=vpos(m, &(incidence[start]), cardinal[n]);
      if(row==-1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
      cefmo.G1.packed[pointer[n]+row]-= factor*P_g *sumx;
      cefmo.G2.packed[pointer[n]+row]-= factor*P_g *sumy;
      row=vpos(mm, &(incidence[start]), cardinal[n]);
      if(row==-1) check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
      cefmo.G1.packed[pointer[n]+row]+= factor*P_g *sumx;
      cefmo.G2.packed[pointer[n]+row]+= factor*P_g *sumy;
//      }
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpGradient2D_CQN1xCQP0(const mesh_t & mesh, cefmo_t & cefmo, parameter_t & data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
 edge x mid-point version

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

{
  int status;

  status=SpGradient2D_CQN1xCQP0_02(mesh, cefmo, data);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SpDivergence2D_CQN1xCQP0(mesh_t & mesh, cefmo_t & cefmo, parameter_t & data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/

{
  int i, j, k, l, m, n, col;
  int status;
  int udim,zdim;
  double S, Lmm, Cd, h;
  double sumx,sumy;
//  double C,factor,pseudo;
  ordering_t ordering;
  int *pointer,*incidence,*cardinal;

  quadrangle_t quadrangle;

  discretisation_t & u_descriptor=*cefmo.u_descriptor;
  discretisation_t & z_descriptor=*cefmo.z_descriptor;

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);
  
  for(m = 0; m < mesh.nquadrangles; m++) {
    quadrangle = mesh.quadrangles[m];
//     pseudo = quadrangle.Area;
//     factor=2.*pseudo/R;
/*-----------------------------------------------------------------------
    compute mass matrix RHS*/
//    C=quadrangle.cosinus;
/**----------------------------------------------------------------------------
    integrale of beta div(Hu) = div(Hu) = sum of Hu.n
-----------------------------------------------------------------------------*/
    for(i = 0; i < z_descriptor.nnpe; i++) {
      l=z_descriptor.NIbE[m][i];
      for(j = 0; j < u_descriptor.nnpe; j++) {
        n=u_descriptor.NIbE[m][j];
        sumx=  quadrangle.l[j]*quadrangle.nx[j];
        sumy=  quadrangle.l[j]*quadrangle.ny[j];
        col=matrix_relative_pos(cefmo.D1.ordering,l,n,1);
        if(col==-1) {
          check_error(-1, "ill-posed matrix", __LINE__, __FILE__, 1);
          }
        size_t pos=cefmo.D1.ordering->pointer[l]+col;
        cefmo.D1.packed[pos]=sumx;
        cefmo.D2.packed[pos]=sumy;
        }
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int WE2D_CQN1xCQP0_SpDivergence(mesh_t & mesh, complex <double> *z, complex <double> *u, complex <double> *v, int nnodes, complex <double> *divergence, const int sign)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  compute non-linear divergence in the continuity equation

  Formulation:

    <div(z u)> = < z u . n>, discontinuity correction applies

  M2xM2 : a cos (wt+p) x b cos (w't+p') = 1/2 ab (cos(wt+p-w't-p') +cos(wt+p+w't+p'))

  cos(p+p')= cos cos - sin sin = real real - imag imag = real(z*z')
  sin(p+p')= cos sin + sin cos = real imag - imag real

  cos(p-p')=cos cos +sin sin= real real + imag imag


  discontinuities:
  
  

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
{
  int i, j, k, l, m, mm, n,status;
  double pseudo;
  complex <double> zz;
  complex <double> sumx,sumy,zint;
  double dzdx,dzdy,C;

  quadrangle_t quadrangle;

  discretisation_t & u_descriptor=mesh.CQN1descriptor;
  discretisation_t & z_descriptor=mesh.CQP0descriptor;
  
  for(m = 0; m < mesh.nquadrangles; m++) {
    quadrangle = mesh.quadrangles[m];
/**----------------------------------------------------------------------------
    integrale of div(Hu) = sum of Hu.n */
    for(i = 0; i < z_descriptor.nnpe; i++) {
      l=z_descriptor.NIbE[m][i];
      for(j = 0; j < u_descriptor.nnpe; j++) {
        n=u_descriptor.NIbE[m][j];
        zz=0;
        for(k=0;k<u_descriptor.nodes[n].nelmts;k++) {
          mm=u_descriptor.nodes[n].elmts[k];
          zz=zz+z[mm];
          }
        zz/= u_descriptor.nodes[n].nelmts;
        sumx=  quadrangle.l[j]*quadrangle.nx[j];
        sumy=  quadrangle.l[j]*quadrangle.ny[j];
        sumx*=  zz*u[n];
        sumy*=  zz*v[n];
/**-----------------------------------------------------------------------
        factor 0.5 from cos a cos b =1/2(cos(a+b) + cos(a-b))*/
        divergence[l] -= 0.5*(sumx+sumy);
        }
      }
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int momentum2d_CQN1xCQP0_SpAdvection(mesh_t & mesh, double *h, complex <double> *u, complex <double> *v,complex <double> *uU, complex <double> *vV, complex <double> *Ax, complex <double> *Ay)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  compute horizontal advection in the conservative momentum equation

    Advection conservative formulation: 

    <beta , div(H u u)> = <div(beta H u u )> - <H u u . grad (beta)>
                        = <div(beta U u )>   - <U u . grad (beta)>

    beta is LGP0, H is LGP1, U is LGP0, u is funny:

    if taken as LGP0:  U u is LGP0, only discontinuities apply

    if taken as ???:  U u is ???


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i, j, k, l, m, n,status;
  triangle_t triangle;
  double pseudo, C;
  complex <double> U,V;
  double sumx,sumy;
  discretisation_t u_descriptor, z_descriptor;

//   u_descriptor=mesh.LGP0descriptor;
//   z_descriptor=mesh.LGP1descriptor;
  check_error(-1, "warning, momentum2d_CQN1xCQP0_SpAdvection not properly implemented yet", __LINE__, __FILE__, 0);

  for(n = 0; n < u_descriptor.nnodes; n++) {
    Ax[n]=0.;
    Ay[n]=0.;
    }
  return(0);

//   for(m = 0; m < mesh.ntriangles; m++) {
//     triangle = mesh.triangles[m];
//     pseudo = triangle.Area;
//     C=triangle.cosinus;
//     U=uU[m]*h[m];
//     V=vV[m]*h[m];
//     for(j = 0; j < 3; j++) {
//       sumx=  triangle.l[j]*triangle.nx[j];
//       sumy=  triangle.l[j]*triangle.ny[j];
// /**-----------------------------------------------------------------------
//       factor 0.5 from discontinuities distribution */
//       Ax[m] += 0.5*sumx*U*u[m]+sumy*U*v[m];
//       Ay[m] += 0.5*sumx*V*u[m]+sumy*V*v[m];
//       }
// /**-----------------------------------------------------------------------
//     factor 0.5 from cos a cos b =1/2(cos(a+b) + cos(a-b))*/
//     Ax[m] *= 0.5;
//     Ay[m] *= 0.5;
// /**-----------------------------------------------------------------------
//     */
//     Ax[m] /= pseudo;
//     Ay[m] /= pseudo;
//     }
// 
//   return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void SpMomentum2D_OBCPatch(mesh_t & mesh, cefmo_t & cefmo, parameter_t & data, tidal_wave wave, complex<double>  *Fu, complex<double>  *Fv)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------

  Spectral momentum solver
  
 -----------------------------------------------------------------------*/
{
  int status;
  int m,n,ntargets;
  int zdim,udim;
  complex<double>  *rhs;
  double theta;
  zmatrix2x2_t A,B,C,Cinv,D,R,Rinv;
  complex<double> nx,ny,tx,ty,unormal;
  cvector2D_t F,S,U;
  complex<double> J=complex<double>(0.,1.);
  complex<double> dzdt,*u,*v;
  vector<paire_t> list;
  vector<paire_t> list2;
  
  discretisation_t *u_descriptor=cefmo.u_descriptor;

  status=paire_dimension(mesh,cefmo.paire,&zdim,&udim);

  for(m=0;m<mesh.nquadrangles;m++) {
    int target, first, count=0;
    for(int i=0;i<u_descriptor->nnpe;i++) {
      n=u_descriptor->NIbE[m][i];
      u_descriptor->nodes[n].code=mesh.edges[n].code;
      if(u_descriptor->nodes[n].code==5) {
        target=i;
        if(count==0) first=i;
        count++;
        }
      }
    if(count==1) {
/**----------------------------------------------------------------------------
      open boundary edge, in an element with only one edge so */
      n=u_descriptor->NIbE[m][target];
      paire_t p(m,target);
      list.push_back(p);
      }
    if(count==2) {
/**----------------------------------------------------------------------------
      open boundary edge, in an element with two edges so */
      n=u_descriptor->NIbE[m][target];
      paire_t p(m,target);
      list2.push_back(p);
      paire_t q(m,first);
      list2.push_back(q);
      }
    }
  
  ntargets=list.size();
  
  u=new complex<double> [ntargets];
  v=new complex<double> [ntargets];
  
  for(int k=0;k<ntargets; k++) {
    paire_t p=list[k];
    m=p.value[0];
    int target=p.value[1];
    n=u_descriptor->NIbE[m][target];
    u[k]=cefmo.state.u[n];
    v[k]=cefmo.state.v[n];
/**----------------------------------------------------------------------------
    necessary to compute continuity equation residuals */
    cefmo.state.u[n]=0;
    cefmo.state.v[n]=0;
    }
  
  rhs =new complex<double> [zdim];
  
  status=matrix_operation(cefmo.D1,cefmo.state.u,rhs,1);
  status=matrix_operation(cefmo.D2,cefmo.state.v,rhs,0);
  
  for(int k=0;k<ntargets; k++) {
    paire_t p=list[k];
    m=p.value[0];
    int target=p.value[1];
    n=u_descriptor->NIbE[m][target];
      
    nx=mesh.quadrangles[m].nx[target];
    ny=mesh.quadrangles[m].ny[target];
    tx=-ny;
    ty= nx;
      
/**----------------------------------------------------------------------------
    continuity equation residual */
    dzdt=cefmo.state.z[m]*J*wave.omega*dph2rps*mesh.CQP0descriptor.massmatrix.packed[m];
    unormal=(-rhs[m]-dzdt)/mesh.edges[n].L;
//    printf("target cell=%d edge=%d %lf %lf\n",m,n,unormal.real(),unormal.imag());
      
    A=zmatrix2x2_t(cefmo.SpMu[n][0],cefmo.SpMu[n][1],cefmo.SpMv[n][0],cefmo.SpMv[n][1]);
      
/**----------------------------------------------------------------------------
    rotation matrix toward tangent,normal components */
    R=zmatrix2x2_t(tx,ty,nx,ny);
/**----------------------------------------------------------------------------
    rotation matrix toward east, north components */
    Rinv=R.inverse();
      
/**----------------------------------------------------------------------------
    momentum equation in tangent,normal components */
    B=R*A*Rinv;
      
/**----------------------------------------------------------------------------
    exchange un and Fn as unknown */
    C=zmatrix2x2_t(B.c[0][0], (complex<double>) 0., B.c[1][0], (complex<double>) -1.);
    Cinv=C.inverse();
      
    F.x=Fu[n]*tx+Fv[n]*ty-B.c[0][1]*unormal;
    F.y=-B.c[1][1]*unormal;
      
    S=Cinv*F;
      
/**----------------------------------------------------------------------------
    S.x is tangent component */
    S.y=unormal;
      
/**----------------------------------------------------------------------------
    retrieve east, north components from tangent,normal components */
    U=Rinv*S;
      
    cefmo.state.u[n]=U.x;
    cefmo.state.v[n]=U.y;
    }

#if 1  
  ntargets=list2.size();
  for(int k=0;k<ntargets; k++) {
    paire_t p=list2[k];
    m=p.value[0];
    int target=p.value[1];
    n=u_descriptor->NIbE[m][target];
    u[k]=cefmo.state.u[n];
    v[k]=cefmo.state.v[n];
/**----------------------------------------------------------------------------
    necessary to compute continuity equation residuals */
    cefmo.state.u[n]=0;
    cefmo.state.v[n]=0;
    }
  
  status=matrix_operation(cefmo.D1,cefmo.state.u,rhs,1);
  status=matrix_operation(cefmo.D2,cefmo.state.v,rhs,0);
  
  for(int k=0;k<ntargets; k++) {
    paire_t p=list2[k];
    m=p.value[0];
    int target=p.value[1];
    n=u_descriptor->NIbE[m][target];
      
    nx=mesh.quadrangles[m].nx[target];
    ny=mesh.quadrangles[m].ny[target];
    tx=-ny;
    ty= nx;
      
/**----------------------------------------------------------------------------
    continuity equation residual */
    dzdt=cefmo.state.z[m]*J*wave.omega*dph2rps*mesh.CQP0descriptor.massmatrix.packed[m];
    unormal=0.5*(-rhs[m]-dzdt)/mesh.edges[n].L;
//    printf("target cell=%d edge=%d\n",m,n);
      
    A=zmatrix2x2_t(cefmo.SpMu[n][0],cefmo.SpMu[n][1],cefmo.SpMv[n][0],cefmo.SpMv[n][1]);
      
/**----------------------------------------------------------------------------
    rotation matrix toward tangent,normal components */
    R=zmatrix2x2_t(tx,ty,nx,ny);
/**----------------------------------------------------------------------------
    rotation matrix toward east, north components */
    Rinv=R.inverse();
      
/**----------------------------------------------------------------------------
    momentum equation in tangent,normal components */
    B=R*A*Rinv;
      
/**----------------------------------------------------------------------------
    exchange un and Fn as unknown */
    C=zmatrix2x2_t(B.c[0][0], (complex<double>) 0., B.c[1][0], (complex<double>) -1.);
    Cinv=C.inverse();
      
    F.x=Fu[n]*tx+Fv[n]*ty-B.c[0][1]*unormal;
    F.y=-B.c[1][1]*unormal;
      
    S=Cinv*F;
      
/**----------------------------------------------------------------------------
    S.x is tangent component */
    S.y=unormal;
      
/**----------------------------------------------------------------------------
    retrieve east, north components from tangent,normal components */
    U=Rinv*S;
      
    cefmo.state.u[n]=U.x;
    cefmo.state.v[n]=U.y;
    }
#endif

#if 0
/**----------------------------------------------------------------------------
  verification */
  status=matrix_operation(cefmo.D1,cefmo.state.u,rhs,1);
  status=matrix_operation(cefmo.D2,cefmo.state.v,rhs,0);
  
  for(int k=0;k<ntargets; k++) {
    paire_t p=list[k];
    m=p.value[0];
    int target=p.value[1];
    n=u_descriptor->NIbE[m][target];
/**----------------------------------------------------------------------------
    continuity equation residual */
    dzdt=cefmo.state.z[m]*J*wave.omega*dph2rps*mesh.CQP0descriptor.massmatrix.packed[m];
    unormal=(-rhs[m]-dzdt)/mesh.edges[n].L;
    printf("target cell=%d edge=%d %lf %lf\n",m,n,unormal.real(),unormal.imag());
    }
    
  for(m=0;m<mesh.nquadrangles;m++) {
    dzdt=cefmo.state.z[m]*J*wave.omega*dph2rps*mesh.CQP0descriptor.massmatrix.packed[m];
    unormal=-rhs[m]-dzdt;
    if(abs(unormal)>1.e-05) {
      printf("target cell=%d  %lf %lf\n",m,unormal.real(),unormal.imag());
      }
    }
#endif

  delete[] u;
  delete[] v;
  delete[] rhs;
  
//  ~list();
  
}

