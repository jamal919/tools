
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Cyril Nguyen       LA/CNRS,    Toulouse, France
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief discretisation projection
*/
/*----------------------------------------------------------------------------*/

#include "tools-structures.h"
#include "fe.h"
#include "zapper.h"
#include "matrix.h"
#include "gauss.h"

static int sproduct2D=0;

//extern int LinearSystem_solve(hypermatrix_t, double *);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void NC2P1_average(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, n, n1, n2;
  edge_t *edges;

  int nelt, nedges, nndes;
  double *w;

  nndes = mesh.nvtxs;
  nelt = mesh.ntriangles;
  edges = mesh.edges;
  nedges = mesh.nedges;

  w = new double[nndes];
  for(n = 0; n < nndes; n++)
    w[n] = 0;
  for(n = 0; n < nndes; n++)
    out[n] = 0;

  for(i = 0; i < nedges; i++) {
    n1 = edges[i].extremity[0];
    n2 = edges[i].extremity[1];
    out[n1] += in[i];
    w[n1]++;
    out[n2] += in[i];
    w[n2]++;
  }

  for(n = 0; n < nndes; n++)
    out[n] /= w[n];
  zaparr(w);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void NC2P1_NQv1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, l, n;
  edge_t *edges;

  int nelts, nedges, nndes;
  double *w;
  double p[3], q[3];
  triangle_t elt;

  nndes = mesh.nvtxs;
  nelts = mesh.ntriangles;
  edges = mesh.edges;
  nedges = mesh.nedges;

  w = new double[nndes];
  for(n = 0; n < nndes; n++)
    w[n] = 0;
  for(n = 0; n < nndes; n++)
    out[n] = 0;

  for(l = 0; l < nelts; l++) {
    elt = mesh.triangles[l];
    for(i = 0; i < 3; i++) {
      n = elt.edges[i];
      q[i] = in[n];
    }
    for(i = 0; i < 3; i++) {
      n = elt.vertex[i];
      for(j = 0; j < 3; j++)
        p[j] = 0;
      p[i] = 1.;
      out[n] += 2. * elt.TrueArea * fe_integraleLGP1xNCP1_2D(p, q);
      w[n] += elt.TrueArea / 3.;
    }
  }

  for(n = 0; n < nndes; n++)
    out[n] /= w[n];
  zaparr(w);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int NC2P1_NQv2(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, l, n;
  edge_t *edges;

  int nelts, nedges, nndes;
  double *w;
  T p[3], q[3];
  triangle_t elt;

  nndes = mesh.nvtxs;
  nelts = mesh.ntriangles;
  edges = mesh.edges;
  nedges = mesh.nedges;

  w = new double[nndes];
  for(n = 0; n < nndes; n++)
    w[n] = 0;
  for(n = 0; n < nndes; n++)
    out[n] = 0;

  for(l = 0; l < nelts; l++) {
    elt = mesh.triangles[l];
    for(i = 0; i < 3; i++) {
      n = elt.edges[i];
      q[i] = in[n];
      }
    for(i = 0; i < 3; i++) {
      p[i] = -q[i];
      for(j = 0; j < 3; j++)
        if(j != i)
          p[i] +=  q[j];
      n = elt.vertex[i];
      out[n] +=  (T) elt.TrueArea * p[i] / (T) 3.;
      w[n]   +=  elt.TrueArea / 3.;
      }
    }

  for(n = 0; n < nndes; n++)
    out[n] /= (T) w[n];
  zaparr(w);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int NC2P1_EI(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, l, n, status;

  double p[3], q[3];
  triangle_t elt;
  int neq, job;

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = 0;

  for(l = 0; l < mesh.ntriangles; l++) {
    elt = mesh.triangles[l];
    for(i = 0; i < 3; i++) {
      n = elt.edges[i];
      q[i] = in[n];
      }
    for(i = 0; i < 3; i++) {
      n = elt.vertex[i];
      for(j = 0; j < 3; j++)
        p[j] = 0;
      p[i] = 1.;
      out[n] += 2. * elt.TrueArea * fe_integraleLGP1xNCP1_2D(p, q);
      }
    }

  neq = mesh.nvtxs;
  job = 0;

  if(mesh.LGP1descriptor.massmatrix.packed==0) {
    check_error(-1, "LGP1descriptor.massmatrix not initialized", __LINE__, __FILE__, 1);
    }
  status=LinearSystem_solve(mesh.LGP1descriptor.massmatrix, out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int NC2P1_EI(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double *a=0,*b=0;
  int n, status;
  size_t size_a,size_b;
  
  size_a=mesh.NCP1descriptor.nnodes;
  size_b=mesh.LGP1descriptor.nnodes;
  
  poc_copy(a,in,size_a);
  b=new double[size_b];
  
  status=NC2P1_EI(a, b, mesh);
  
  for (n=0;n< size_b;n++) out[n]=b[n];
  
  delete[] a;
  delete[] b;
  
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int NC2P1_EI_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double *a,*b,*c;
  int n, status;
  size_t size_a,size_b;
  
  size_a=mesh.NCP1descriptor.nnodes;
  size_b=mesh.LGP1descriptor.nnodes;
  
  a=new double[size_a];
  b=new double[size_b];
  c=new double[size_b];
  
  for (n=0;n< size_a;n++) a[n]=in[n].real();
  status=NC2P1_EI(a, b, mesh);
  for (n=0;n< size_a;n++) a[n]=in[n].imag();
  status=NC2P1_EI(a, c, mesh);
  
  for (n=0;n< size_b;n++) out[n]=complex< T >(b[n],c[n]);
  
  delete[] a;
  delete[] b;
  delete[] c;
  
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int NC2P1_EI(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  
  status=NC2P1_EI_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int NC2P1_EI(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  
  status=NC2P1_EI_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int projection_NCP1xLGP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------

  Project NCP1 field on LGP1 field

-----------------------------------------------------------------------------*/
{
  int status;
  
  if(in==0) return(-1);

  switch (sproduct2D) {
    case (NQUAD):
/*------------------------------------------------------------------------------
      Nodal quadrature */
      status=NC2P1_NQv2(in, out, mesh);
      break;

    case (INTGL):
/*------------------------------------------------------------------------------
      True integral */
      status=NC2P1_EI(in, out, mesh);
      break;
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int projection_NCP1xLGP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=projection_NCP1xLGP1_template(in,out,mesh);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int projection_NCP1xLGP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=projection_NCP1xLGP1_template(in,out,mesh);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int projection_NCP1xLGP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=projection_NCP1xLGP1_template(in,out,mesh);
  
  return result;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int projection_NCP1xLGP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=projection_NCP1xLGP1_template(in,out,mesh);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int projection_NCP1xLGP0_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Project NCP1 field on LGP0 field

  uses integrale norm

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, m, n;
  T q[3];
  triangle_t elt;

  if(in==0) return(-1);

  for(m = 0; m < mesh.ntriangles; m++)
    out[m] = 0;

  for(m = 0; m < mesh.ntriangles; m++) {
    elt = mesh.triangles[m];
    for(i = 0; i < 3; i++) {
      n = elt.edges[i];
      q[i] = in[n];
      }
/*-----------------------------------------------------------------------------
    do not multiply by area as it would be divided later by the same quantity*/
    out[m] += (T) 2. * fe_integraleNCP1_2D(q);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1xLGP0(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------

  Project NCP1 field on LGP0 field

  uses integrale norm

-----------------------------------------------------------------------------*/
{
  int m;

  if(in==0) return(-1);

  for(m = 0; m < mesh.ntriangles; m++)
    out[m] = 0;

#pragma omp parallel for //if(gOPENMP_nCPUs>1)
  for(m = 0; m < mesh.ntriangles; m++) {
    int i,n;
    double q[3];
    triangle_t elt = mesh.triangles[m];
    for(i = 0; i < 3; i++) {
      n = elt.edges[i];
      q[i] = in[n];
      }
/**----------------------------------------------------------------------------
    do not multiply by area as it would be divided later by the same quantity*/
    out[m] += 2. * fe_integraleNCP1_2D(q);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP1xLGP0_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n;
  T sum,p[3];

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      p[k]=in[n];
      }
    sum=0;
    for(k=0;k<3;k++) {
      sum+=p[k];
      }
    out[m]=sum/(T) 3.0;
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int projection_LGP1xLGP0(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xLGP0_template(in,out,mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int projection_LGP1xLGP0(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xLGP0_template(in,out,mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int projection_LGP1xLGP0(complex<double>  *in, complex<double>  *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xLGP0_template(in,out,mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xLGP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double precision vector projection from LGP0 nodes onto LGP1 nodes */
{
  int k,m,n,status;
  double pseudo,c[3];
  extern int print_vector(char* file, mesh_t mesh, int discretisation, double *val);

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = 0;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      c[k]=mesh.vertices[n].c;
      }
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      out[n] += 2.0 * pseudo * fe_integraleLGP1xLGP1_2D(c,k) * in[m];
      }
    }
   
  status= LinearSystem_solve(mesh.LGP1descriptor.massmatrix, out);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xLGP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double *a=0,*b=0;
  int n, status;
  size_t size_a,size_b;
  
  size_a=mesh.LGP0descriptor.nnodes;
  if(size_a==0) return(-1);
  size_b=mesh.LGP1descriptor.nnodes;
  if(size_b==0) return(-1);
  
  poc_copy(a,in,size_a);
  b=new double[size_b];
  
  status=projection_LGP0xLGP1(a, b, mesh);
  
  for (n=0;n< size_b;n++) out[n]=b[n];
  
  delete[] a;
  delete[] b;
  
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xLGP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes onto LGP1 nodes */
{
  int k,m,n,status;
  double pseudo,C,c[3],sum;
  complex<double> zz;
  double *real,*imag;
//  extern int LinearSystem_solve(hypermatrix_t &, double *);

  discretisation_t descriptor=mesh.LGP1descriptor;

  if(in==0) return(-1);

  real=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    real[n] = 0;
  imag=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
      c[k]=mesh.vertices[n].c;
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      sum=fe_integraleLGP1xLGP1_2D(c,k);
      if(isnan(in[m].real())) {
        printf("projection_LGP0xLGP1 : troubles...\n");
        }
/**----------------------------------------------------------------------------
        bug? cosinus is already accounted for in integrale*/
//      zz = 2.0 * pseudo * C * sum * in[m];
      zz = 2.0 * pseudo * sum * in[m];
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,real);
  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,imag);

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = complex<double> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_LGP0xLGP1_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes onto LGP1 nodes */
{
  int k,m,n,status;
  double pseudo,C,c[3],sum;
  complex<double> zz;
  double *real,*imag;

  discretisation_t descriptor=mesh.LGP1descriptor;

  if(in==0) return(-1);

  real=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    real[n] = 0;
  imag=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
      c[k]=mesh.vertices[n].c;
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      sum=fe_integraleLGP1xLGP1_2D(c,k);
      if(isnan(in[m].real())) {
        printf("projection_LGP0xLGP1 : troubles...\n");
        }
/**----------------------------------------------------------------------------
        bug? cosinus is already accounted for in integrale*/
//      zz = 2.0 * pseudo * C * sum * in[m];
      zz = (T) (2.0 * pseudo * sum) * in[m];
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,real);
  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,imag);

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = complex<T> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xLGP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes onto LGP1 nodes */
{
  int status;
  
  status=projection_LGP0xLGP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_LGP0xNCP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to NCP1 nodes */
{
  int m,status;
  int  m1,m2;
  T *weight;

  if(in==0) return(-1);

  for(size_t n = 0; n < mesh.nedges; n++)
    out[n] = 0;

  weight=new T[mesh.nedges];

  for(size_t n = 0; n < mesh.nedges; n++)
    weight[n] = 0;

//#pragma omp parallel for if(gOPENMP_nCPUs>1)
  for (m=0;m<mesh.ntriangles;m++) {
    int k,n;
    double pseudo=mesh.triangles[m].Area;
    double c[3],w;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      c[k]=mesh.vertices[n].c;
      }
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].edges[k];
      w=pseudo*fe_integraleLGP1xNCP1_2D(c,k);
//       out[n]    += 2.0 * pseudo * fe_integraleLGP1xNCP1_2D(c,k) * in[m];
//       weight[n] += 2.0 * pseudo * fe_integraleLGP1xNCP1_2D(c,k);
      out[n]    += (T) w * in[m];
      weight[n] += (T) w;
      }
    }

  for(size_t n = 0; n < mesh.nedges; n++)
    out[n]/=weight[n];

  delete[] weight;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xNCP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to NCP1 nodes */
{
  int status;
  
  status=projection_LGP0xNCP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xNCP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to NCP1 nodes */
{
  int status;
  
  status=projection_LGP0xNCP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xNCP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform float complex vector projection from LGP0 nodes to NCP1 nodes */
{
  int status;
  
  status=projection_LGP0xNCP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xNCP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to NCP1 nodes */
{
  int status;
  
  status=projection_LGP0xNCP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP1xDGP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,nLGP1,status=0;

  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<mesh.DGP1descriptor.nnpe;k++) {
      n=mesh.DGP1descriptor.NIbE[m][k];
      nLGP1=mesh.triangles[m].vertex[k];
      out[n]   = in[nLGP1];
      }
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xDGP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xDGP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xDGP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xDGP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xDGP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xDGP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xDGP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xDGP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP1xNCP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,n1,n2,status;

  if(in==0) return(-1);

  for (n=0;n<mesh.nedges;n++) {
    n1=mesh.edges[n].extremity[0];
    n2=mesh.edges[n].extremity[1];
    out[n] = (T) 0.5*(in[n1]+in[n2]);
    }

  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xNCP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xNCP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xNCP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xNCP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xNCP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xNCP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xNCP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xNCP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP1xDNP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,n1,n2,status;

  if(in==0) return(-1);

  for (n=0;n<mesh.DNP1descriptor.nnodes;n++) {
    n1=mesh.DNP1descriptor.nodes[n].vtces[0];
    n2=mesh.DNP1descriptor.nodes[n].vtces[1];
    out[n] = (T) 0.5*(in[n1]+in[n2]);
    }

  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xDNP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xDNP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xDNP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xDNP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xDNP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xDNP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xDNP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP1xDNP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xLGP2(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double real vector projection from LGP0 nodes to LGP2 nodes */
{
  int k,m,n,status;
  double pseudo/*,c[3]*/,C;
  double *buffer;
  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  buffer=new double[descriptor.nnodes];

  for(n = 0; n < descriptor.nnodes; n++)
    buffer[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    C=mesh.triangles[m].cosinus;
    for(k=0;k<6;k++) {
      n=descriptor.NIbE[m][k];
      buffer[n] += 2.0 * pseudo * C * fe_LGP2_2D[k] * in[m];
      }
    }

  status = LinearSystem_solve(descriptor.massmatrix,buffer);
  for(n = 0; n < descriptor.nnodes; n++)
    out[n] = (float) buffer[n];

  delete[] buffer;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xLGP2(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double real vector projection from LGP0 nodes to LGP2 nodes */
{
  int k,m,n,status;
  double pseudo/*,c[3]*/,C;
  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  for(n = 0; n < descriptor.nnodes; n++)
    out[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
//     for(k=0;k<3;k++) {
//       n=mesh.triangles[m].vertex[k];
//       c[k]=mesh.vertices[n].c;
//       }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<6;k++) {
      n=descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_LGP2_2D[k] * in[m];
      }
    }

  status = LinearSystem_solve(descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xLGP2(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to LGP2 nodes */
{
  int n,status;
  double *in_real,*in_imag;
  double *out_real,*out_imag;
  discretisation_t in_descriptor =mesh.LGP0descriptor;
  discretisation_t out_descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  in_real=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_real[n] = in[n].real();
  in_imag=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_imag[n] = in[n].imag();

  out_real=new double[out_descriptor.nnodes];
  out_imag=new double[out_descriptor.nnodes];

  status= projection_LGP0xLGP2(in_real, out_real, mesh);
  status= projection_LGP0xLGP2(in_imag, out_imag, mesh);

  for(n = 0; n < out_descriptor.nnodes; n++)
    out[n] = complex<double> (out_real[n],out_imag[n]);

  delete[] in_real;
  delete[] in_imag;
  delete[] out_real;
  delete[] out_imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_LGP0xLGP2_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to LGP2 nodes */
{
  int n,status;
  double *in_real,*in_imag;
  double *out_real,*out_imag;
  discretisation_t in_descriptor =mesh.LGP0descriptor;
  discretisation_t out_descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  in_real=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_real[n] = in[n].real();
  in_imag=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_imag[n] = in[n].imag();

  out_real=new double[out_descriptor.nnodes];
  out_imag=new double[out_descriptor.nnodes];

  status= projection_LGP0xLGP2(in_real, out_real, mesh);
  status= projection_LGP0xLGP2(in_imag, out_imag, mesh);

  for(n = 0; n < out_descriptor.nnodes; n++)
    out[n] = complex<T> (out_real[n],out_imag[n]);

  delete[] in_real;
  delete[] in_imag;
  delete[] out_real;
  delete[] out_imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xLGP2(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to LGP2 nodes */
{
  int n,status;
  
  status=projection_LGP0xLGP2_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP0xDGP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double real vector projection from LGP0 nodes to DGP1 nodes */
{
  int k,m,n,status=0;
  double pseudo,C;
  discretisation_t descriptor=mesh.DGP1descriptor;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=descriptor.NIbE[m][k];
      out[n] = in[m];
      }
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xDGP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double real vector projection from LGP0 nodes to DGP1 nodes */
{
  int status;

  status =projection_LGP0xDGP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xDGP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double real vector projection from LGP0 nodes to DGP1 nodes */
{
  int status;

  status =projection_LGP0xDGP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xDGP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to LGP2 nodes */
{
  int status;

  status =projection_LGP0xDGP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xDGP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to LGP2 nodes */
{
  int status;

  status =projection_LGP0xDGP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP0xDNP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double real vector projection from LGP0 nodes to DGP1 nodes */
{
  int k,m,n,status;
  discretisation_t descriptor=mesh.DNP1descriptor;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    for(k=0;k<3;k++) {
      n=descriptor.NIbE[m][k];
      out[n] = in[m];
      }
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xDNP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double real vector projection from LGP0 nodes to DGP1 nodes */
{
  int status;

  status =projection_LGP0xDNP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xDNP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double real vector projection from LGP0 nodes to DGP1 nodes */
{
  int status;

  status =projection_LGP0xDNP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0xDNP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to LGP2 nodes */
{
  int status;

  status =projection_LGP0xDNP1_template(in, out, mesh);
  
  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP2xLGP0_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n;
  T sum;
  double w[6],weight;

  if(in==0) return(-1);

  weight=0;
  for(k=0;k<6;k++) {
    w[k]=2.0*fe_LGP2_2D[k];
    weight+=w[k];
    }

  for (m=0;m<mesh.ntriangles;m++) {
    sum=0;
    for(k=1;k<6;k+=2) {
      n=mesh.LGP2descriptor.NIbE[m][k];
      sum+=in[n];
      }
    out[m]=sum/(T) 3.0;
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xLGP0(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n;
  double sum/*,p[6]*/,w[6],weight;

  if(in==0) return(-1);

  weight=0;
  for(k=0;k<6;k++) {
    w[k]=2.0*fe_LGP2_2D[k];
    weight+=w[k];
    }

  for (m=0;m<mesh.ntriangles;m++) {
//     for(k=0;k<6;k++) {
//       n=mesh.LGP2descriptor.NIbE[m][k];
//       p[k]=in[n];
//       }
    sum=0;
//     for(k=0;k<6;k++) {
//       sum+=w[k]*p[k];
//       }
//     out[m]=sum/weight;
    for(k=1;k<6;k+=2) {
      n=mesh.LGP2descriptor.NIbE[m][k];
      sum+=in[n];
      }
    out[m]=sum/3.0;
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xLGP2(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[3];
  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  for(n = 0; n < descriptor.nnodes; n++)
    out[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(k,p);
      }
    }
  status = LinearSystem_solve(descriptor.massmatrix,out);
  return(status);
}

//*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xLGP2(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[3];
  double *buffer;

  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  buffer=new double[descriptor.nnodes];

  for(n = 0; n < descriptor.nnodes; n++)
    buffer[n] = 0;


  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(k,p);
      }
    }
  status = LinearSystem_solve(descriptor.massmatrix,buffer);
  for(n = 0; n < descriptor.nnodes; n++)
    out[n] = (float) buffer[n];

  delete[] buffer;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xLGP2(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP1 nodes to LGP2 nodes */
{
  int n,status;
  double *in_real,*in_imag;
  double *out_real,*out_imag;
  discretisation_t in_descriptor =mesh.LGP1descriptor;
  discretisation_t out_descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  in_real=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_real[n] = in[n].real();
  in_imag=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_imag[n] = in[n].imag();

  out_real=new double[out_descriptor.nnodes];
  out_imag=new double[out_descriptor.nnodes];

  status= projection_LGP1xLGP2(in_real, out_real, mesh);
  status= projection_LGP1xLGP2(in_imag, out_imag, mesh);

  for(n = 0; n < out_descriptor.nnodes; n++)
    out[n] = complex<double> (out_real[n],out_imag[n]);

  delete[] in_real;
  delete[] in_imag;
  delete[] out_real;
  delete[] out_imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_LGP1xLGP2_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP1 nodes to LGP2 nodes */
{
  int n,status;
  double *in_real,*in_imag;
  double *out_real,*out_imag;
  discretisation_t in_descriptor =mesh.LGP1descriptor;
  discretisation_t out_descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  in_real=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_real[n] = in[n].real();
  in_imag=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_imag[n] = in[n].imag();

  out_real=new double[out_descriptor.nnodes];
  out_imag=new double[out_descriptor.nnodes];

  status= projection_LGP1xLGP2(in_real, out_real, mesh);
  status= projection_LGP1xLGP2(in_imag, out_imag, mesh);

  for(n = 0; n < out_descriptor.nnodes; n++)
    out[n] = complex<T> (out_real[n],out_imag[n]);

  delete[] in_real;
  delete[] in_imag;
  delete[] out_real;
  delete[] out_imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1xLGP2(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP1 nodes to LGP2 nodes */
{
  int status;
  
  status=projection_LGP1xLGP2_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_NCP1xDNP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,n1,n2,status;
  discretisation_t descriptor=mesh.DNP1descriptor;

  if(in==0) return(-1);

  for(n = 0; n < mesh.nedges; n++) {
    switch(mesh.edges[n].nshared) {
      case 1:
        i=vpos(n,mesh.triangles[mesh.edges[n].shared[0]].edges,3);
        n1=descriptor.NIbE[mesh.edges[n].shared[0]][i];
        out[n1]=in[n];
        break;;
      case 2:
        i=vpos(n,mesh.triangles[mesh.edges[n].shared[0]].edges,3);
        j=vpos(n,mesh.triangles[mesh.edges[n].shared[1]].edges,3);
        n1=descriptor.NIbE[mesh.edges[n].shared[0]][i];
        n2=descriptor.NIbE[mesh.edges[n].shared[1]][j];
        out[n1]=in[n];
        out[n2]=in[n];
        break;;
      default:
       check_error(-1, "anomaly in edges structure", __LINE__, __FILE__, 1);
       break;
      }
    }
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1xDNP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=projection_NCP1xDNP1_template(in, out,  mesh);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1xDNP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=projection_NCP1xDNP1_template(in, out,  mesh);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1xDNP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=projection_NCP1xDNP1_template(in, out,  mesh);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1xDNP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=projection_NCP1xDNP1_template(in, out,  mesh);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//  int projection_NCP1xLGP2(double *in, double *out, mesh_t mesh, discretisation_t descriptor)
  int projection_NCP1xLGP2(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[3];
  double *buffer;

  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  buffer=new double[descriptor.nnodes];

  for(n = 0; n < descriptor.nnodes; n++)
    buffer[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].edges[k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(k,p);
      }
    }
  status = LinearSystem_solve(descriptor.massmatrix,buffer);
  for(n = 0; n < descriptor.nnodes; n++)
    out[n] = buffer[n];

  delete[] buffer;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//  int projection_NCP1xLGP2(double *in, double *out, mesh_t mesh, discretisation_t descriptor)
  int projection_NCP1xLGP2(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[3];

  discretisation_t descriptor=mesh.LGP2descriptor;

  for(n = 0; n < descriptor.nnodes; n++)
    out[n] = 0;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].edges[k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(k,p);
      }
    }
  status = LinearSystem_solve(descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1xLGP2(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to LGP2 nodes */
{
  int n,status;
  double *in_real,*in_imag;
  double *out_real,*out_imag;
  discretisation_t in_descriptor =mesh.NCP1descriptor;
  discretisation_t out_descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  in_real=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_real[n] = in[n].real();
  in_imag=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_imag[n] = in[n].imag();

  out_real=new double[out_descriptor.nnodes];
  out_imag=new double[out_descriptor.nnodes];

  status= projection_NCP1xLGP2(in_real, out_real, mesh);
  status= projection_NCP1xLGP2(in_imag, out_imag, mesh);

  for(n = 0; n < out_descriptor.nnodes; n++)
    out[n] = complex<double> (out_real[n],out_imag[n]);

  delete[] in_real;
  delete[] in_imag;
  delete[] out_real;
  delete[] out_imag;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_NCP1xLGP2_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Perform double complex vector projection from LGP0 nodes to LGP2 nodes */
{
  int n,status;
  double *in_real,*in_imag;
  double *out_real,*out_imag;
  discretisation_t in_descriptor =mesh.NCP1descriptor;
  discretisation_t out_descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  in_real=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_real[n] = in[n].real();
  in_imag=new double[in_descriptor.nnodes];
  for(n = 0; n < in_descriptor.nnodes; n++)
    in_imag[n] = in[n].imag();

  out_real=new double[out_descriptor.nnodes];
  out_imag=new double[out_descriptor.nnodes];

  status= projection_NCP1xLGP2(in_real, out_real, mesh);
  status= projection_NCP1xLGP2(in_imag, out_imag, mesh);

  for(n = 0; n < out_descriptor.nnodes; n++)
    out[n] = complex<T> (out_real[n],out_imag[n]);

  delete[] in_real;
  delete[] in_imag;
  delete[] out_real;
  delete[] out_imag;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1xLGP2(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_NCP1xLGP2_template(in,out,mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1xLGP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[6];

  discretisation_t descriptor=mesh.DGP1descriptor;

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = 0;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP1xLGP1_2D(p,k);
      }
    }
  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1xLGP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double *a=0,*b=0;
  int n, status;
  size_t size_a,size_b;
  
  size_a=mesh.DGP1descriptor.nnodes;
  size_b=mesh.LGP1descriptor.nnodes;
  
  poc_copy(a,in,size_a);
  b=new double[size_b];
  
  status=projection_DGP1xLGP1(a, b, mesh);
  
  for (n=0;n< size_b;n++) out[n]=b[n];
  
  delete[] a;
  delete[] b;
  
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_DGP1xLGP1_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[3];
  double *real,*imag;

  discretisation_t descriptor=mesh.DGP1descriptor;

  if(in==0) return(-1);

  real=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    real[n] = 0;
  imag=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      zz = 2.0 * pseudo * C * fe_integraleLGP1xLGP1_2D(p,k);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,real);
  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,imag);

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = complex<T> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1xLGP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_DGP1xLGP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1xLGP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_DGP1xLGP1_template(in, out, mesh);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1xLGP2(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  double zz,p[3];
  double *tmp;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.DGP1descriptor;

  if(in==0) return(-1);

  tmp=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    tmp[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(k,p);
      tmp[n]+=zz;
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,tmp);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = (float) tmp[n];

  delete[] tmp;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1xLGP2(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[3];

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.DGP1descriptor;

  if(in==0) return(-1);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(k,p);
      }
    }
  status = LinearSystem_solve(dst_descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1xLGP2(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[3];
  double *real,*imag;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.DGP1descriptor;

  if(in==0) return(-1);

  real=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(k,p);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,real);
  status = LinearSystem_solve(dst_descriptor.massmatrix,imag);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = complex<double> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_DGP1xLGP2_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[3];
  double *real,*imag;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.DGP1descriptor;

  if(in==0) return(-1);

  real=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(k,p);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,real);
  status = LinearSystem_solve(dst_descriptor.massmatrix,imag);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = complex<T> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1xLGP2(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=projection_DGP1xLGP2(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xLGP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[6];

  discretisation_t descriptor=mesh.DNP1descriptor;

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = 0;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP1xNCP1_2D(k,p);
      }
    }
  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xLGP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double *a=0,*b=0;
  int n, status;
  size_t size_a,size_b;
  
  size_a=mesh.DNP1descriptor.nnodes;
  size_b=mesh.LGP1descriptor.nnodes;
  
  poc_copy(a,in,size_a);
  b=new double[size_b];
  
  status=projection_DNP1xLGP1(a, b, mesh);
  
  for (n=0;n< size_b;n++) out[n]=b[n];
  
  delete[] a;
  delete[] b;
  
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template < typename T >  int projection_DNP1xLGP1_template(complex<T> *in, complex<T> *out, mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[3];
  double *real,*imag;

  discretisation_t descriptor=mesh.DNP1descriptor;

  if(in==0) return(-1);

  real=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    real[n] = 0;
  imag=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      zz = 2.0 * pseudo * C * fe_integraleLGP1xNCP1_2D(k,p);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,real);
  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,imag);

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = complex<T> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xLGP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_DNP1xLGP1_template(in, out, mesh);
  
  return(status);
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xLGP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_DNP1xLGP1_template(in, out, mesh);
  
  return(status);
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_DNP1xNCP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,n1,n2,status;
  discretisation_t descriptor=mesh.DNP1descriptor;
  T mean;

  if(in==0) return(-1);

  for(n = 0; n < mesh.nedges; n++) {
    switch(mesh.edges[n].nshared) {
      case 1:
        i=vpos(n,mesh.triangles[mesh.edges[n].shared[0]].edges,3);
        n1=descriptor.NIbE[mesh.edges[n].shared[0]][i];
        out[n]=in[n1];
        break;
      case 2:
        i=vpos(n,mesh.triangles[mesh.edges[n].shared[0]].edges,3);
        j=vpos(n,mesh.triangles[mesh.edges[n].shared[1]].edges,3);
        n1=descriptor.NIbE[mesh.edges[n].shared[0]][i];
        n2=descriptor.NIbE[mesh.edges[n].shared[1]][j];
        mean=0.5*(in[n1]+in[n2]);
        out[n]=mean;
        break;
      default:
       check_error(-1, "anomaly in edges structure", __LINE__, __FILE__, 1);
       break;
      }
    }

  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xNCP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=projection_DNP1xNCP1_template(in, out,  mesh);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xNCP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

//  status=projection_DNP1xNCP1_template(in, out,  mesh);
   check_error(-1, "projection not yet implemented", __LINE__, __FILE__, 1);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xNCP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

//  status=projection_DNP1xNCP1_template(in, out,  mesh);
   check_error(-1, "projection not yet implemented", __LINE__, __FILE__, 1);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xNCP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

//  status=projection_DNP1xNCP1_template(in, out,  mesh);
   check_error(-1, "projection not yet implemented", __LINE__, __FILE__, 1);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xLGP2(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  double zz,p[3];
  double *tmp;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.DNP1descriptor;

  if(in==0) return(-1);

  tmp=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    tmp[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(k,p);
      tmp[n]+=zz;
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,tmp);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = (float) tmp[n];

  delete[] tmp;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xLGP2(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[3];

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.DNP1descriptor;

  if(in==0) return(-1);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(k,p);
      }
    }
  status = LinearSystem_solve(dst_descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xLGP2(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[3];
  double *real,*imag;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.DNP1descriptor;

  if(in==0) return(-1);

  real=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(k,p);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,real);
  status = LinearSystem_solve(dst_descriptor.massmatrix,imag);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = complex<double> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_DNP1xLGP2_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[3];
  double *real,*imag;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.DNP1descriptor;

  if(in==0) return(-1);

  real=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(k,p);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,real);
  status = LinearSystem_solve(dst_descriptor.massmatrix,imag);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = complex<T> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1xLGP2(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_DNP1xLGP2_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xLGP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[6];

  discretisation_t descriptor=mesh.LGP2descriptor;

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = 0;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(p,k);
      }
    }
  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xLGP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double *a=0,*b=0;
  int n, status;
  size_t size_a,size_b;
  
  size_a=mesh.LGP2descriptor.nnodes;
  size_b=mesh.LGP1descriptor.nnodes;
  
  poc_copy(a,in,size_a);
  b=new double[size_b];
  
  status=projection_LGP2xLGP1(a, b, mesh);
  
  for (n=0;n< size_b;n++) out[n]=b[n];
  
  delete[] a;
  delete[] b;
  
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_LGP2xLGP1_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[6];
  double *real,*imag;
  extern complex<double> fe_integraleLGP2xLGP1_2D(complex<double> *p, int j);

  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  real=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    real[n] = 0;
  imag=new double[mesh.nvtxs];
  for(n = 0; n < mesh.nvtxs; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.triangles[m].vertex[k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(p,k);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,real);
  status = LinearSystem_solve(mesh.LGP1descriptor.massmatrix,imag);

  for(n = 0; n < mesh.nvtxs; n++)
    out[n] = complex<T> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xLGP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP2xLGP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xLGP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP2xLGP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xDGP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[6];
  discretisation_t descriptor=mesh.LGP2descriptor;
  double *buffer;

  if(in==0) return(-1);

  buffer=new double[mesh.DGP1descriptor.nnodes];

  for(n = 0; n < mesh.DGP1descriptor.nnodes; n++)
    buffer[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.DGP1descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(p,k);
      }
    }
  status = LinearSystem_solve(mesh.DGP1descriptor.massmatrix,buffer);
  for(n = 0; n < mesh.DGP1descriptor.nnodes; n++)
    out[n] = buffer[n];

  delete[] buffer;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xDGP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[6];
  discretisation_t src_descriptor=mesh.LGP2descriptor;
  discretisation_t dst_descriptor=mesh.DGP1descriptor;

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = 0;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(p,k);
      }
    }
  status = LinearSystem_solve(dst_descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xDGP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[6];
  double *real,*imag;
  extern complex<double> fe_integraleLGP2xLGP1_2D(complex<double> *p, int j);

  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  real=new double[mesh.DGP1descriptor.nnodes];
  for(n = 0; n < mesh.DGP1descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[mesh.DGP1descriptor.nnodes];
  for(n = 0; n < mesh.DGP1descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.DGP1descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(p,k);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(mesh.DGP1descriptor.massmatrix,real);
  status = LinearSystem_solve(mesh.DGP1descriptor.massmatrix,imag);

  for(n = 0; n < mesh.DGP1descriptor.nnodes; n++)
    out[n] = complex<double> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_LGP2xDGP1_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[6];
  double *real,*imag;
  extern complex<double> fe_integraleLGP2xLGP1_2D(complex<double> *p, int j);

  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  real=new double[mesh.DGP1descriptor.nnodes];
  for(n = 0; n < mesh.DGP1descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[mesh.DGP1descriptor.nnodes];
  for(n = 0; n < mesh.DGP1descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.DGP1descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xLGP1_2D(p,k);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(mesh.DGP1descriptor.massmatrix,real);
  status = LinearSystem_solve(mesh.DGP1descriptor.massmatrix,imag);

  for(n = 0; n < mesh.DGP1descriptor.nnodes; n++)
    out[n] = complex<double> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xDGP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP2xDGP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xDNP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[6];
  discretisation_t descriptor=mesh.LGP2descriptor;
  double *buffer;

  if(in==0) return(-1);

  buffer=new double[mesh.DNP1descriptor.nnodes];

  for(n = 0; n < mesh.DNP1descriptor.nnodes; n++)
    buffer[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.DNP1descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(p,k);
      }
    }
  status = LinearSystem_solve(mesh.DNP1descriptor.massmatrix,buffer);
  for(n = 0; n < mesh.DNP1descriptor.nnodes; n++)
    out[n] = buffer[n];

  delete[] buffer;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xDNP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[6];
  discretisation_t src_descriptor=mesh.LGP2descriptor;
  discretisation_t dst_descriptor=mesh.DNP1descriptor;

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = 0;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(p,k);
      }
    }
  status = LinearSystem_solve(dst_descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_LGP2xDNP1_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[6];
  double *real,*imag;
  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  real=new double[mesh.DNP1descriptor.nnodes];
  for(n = 0; n < mesh.DNP1descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[mesh.DNP1descriptor.nnodes];
  for(n = 0; n < mesh.DNP1descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.DNP1descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(p,k);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(mesh.DNP1descriptor.massmatrix,real);
  status = LinearSystem_solve(mesh.DNP1descriptor.massmatrix,imag);

  for(n = 0; n < mesh.DNP1descriptor.nnodes; n++)
    out[n] = complex<double> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xDNP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_LGP2xDNP1_template(in, out, mesh);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xDNP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[6];
  double *real,*imag;
  discretisation_t descriptor=mesh.LGP2descriptor;

  if(in==0) return(-1);

  real=new double[mesh.DNP1descriptor.nnodes];
  for(n = 0; n < mesh.DNP1descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[mesh.DNP1descriptor.nnodes];
  for(n = 0; n < mesh.DNP1descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<descriptor.nnpe;k++) {
      n=descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<3;k++) {
      n=mesh.DNP1descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(p,k);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(mesh.DNP1descriptor.massmatrix,real);
  status = LinearSystem_solve(mesh.DNP1descriptor.massmatrix,imag);

  for(n = 0; n < mesh.DNP1descriptor.nnodes; n++)
    out[n] = complex<double> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_IPGxLGP2(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  double p[16];
  double *tmp;
  gauss_t gauss;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.IPGdescriptor;

  if(in==0) return(-1);
  
  gauss.init(mesh.IPGdescriptor.type);

  tmp=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    tmp[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      tmp[n] = 2.0 * pseudo * C * fe_QIntegraleLGP2xIPG7_2D(k,p,gauss);
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,tmp);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = tmp[n];

  delete[] tmp;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_IPGxLGP2(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  double p[16];
  gauss_t gauss;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.IPGdescriptor;

  if(in==0) return(-1);
  
  gauss.init(mesh.IPGdescriptor.type);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      out[n]+= 2.0 * pseudo * C * fe_QIntegraleLGP2xIPG7_2D(k,p,gauss);
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,out);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_IPGxLGP2(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[16];
  double *real,*imag;
  gauss_t gauss;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.IPGdescriptor;

  if(in==0) return(-1);
  
  gauss.init(mesh.IPGdescriptor.type);

  real=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
//      p[k]=src_descriptor.nodes[n].c;
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_QIntegraleLGP2xIPG7_2D(k,p,gauss);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,real);
  status = LinearSystem_solve(dst_descriptor.massmatrix,imag);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = complex<double> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  
  gauss.destroy();
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_IPGxLGP2_template(complex<T> *in, complex<T> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/;
  complex<double> zz,p[16];
  double *real,*imag;
  gauss_t gauss;

  discretisation_t dst_descriptor=mesh.LGP2descriptor;
  discretisation_t src_descriptor=mesh.IPGdescriptor;

  if(in==0) return(-1);
  
  gauss.init(mesh.IPGdescriptor.type);

  real=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    real[n] = 0;
  imag=new double[dst_descriptor.nnodes];
  for(n = 0; n < dst_descriptor.nnodes; n++)
    imag[n] = 0;

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
//      p[k]=src_descriptor.nodes[n].c;
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      zz = 2.0 * pseudo * C * fe_QIntegraleLGP2xIPG7_2D(k,p,gauss);
      real[n]+=zz.real();
      imag[n]+=zz.imag();
      }
    }

  status = LinearSystem_solve(dst_descriptor.massmatrix,real);
  status = LinearSystem_solve(dst_descriptor.massmatrix,imag);

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = complex<T> (real[n],imag[n]);

  delete[] real;
  delete[] imag;
  
  gauss.destroy();
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_IPGxLGP2(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_IPGxLGP2_template(in, out, mesh);

  return(status);
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP2xNCP1_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status;
  double pseudo,C/*,c[3]*/,p[6];
  discretisation_t src_descriptor=mesh.LGP2descriptor;
  discretisation_t dst_descriptor=mesh.NCP1descriptor;

  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = 0;

  if(in==0) return(-1);

  for (m=0;m<mesh.ntriangles;m++) {
    pseudo=mesh.triangles[m].Area;
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
//      c[k]=mesh.vertices[n].c;
      p[k]=in[n];
      }
    C=mesh.triangles[m].cosinus;
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      out[n] += 2.0 * pseudo * C * fe_integraleLGP2xNCP1_2D(p,k);
      }
    }
  status = LinearSystem_solve(dst_descriptor.massmatrix,out);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xNCP1(double *in, double *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=projection_LGP2xNCP1_template(in, out,  mesh);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xNCP1(float *in, float *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

//  status=projection_DNP1xNCP1_template(in, out,  mesh);
   check_error(-1, "projection not yet implemented", __LINE__, __FILE__, 1);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xNCP1(complex<float> *in, complex<float> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

//  status=projection_DNP1xNCP1_template(in, out,  mesh);
   check_error(-1, "projection not yet implemented", __LINE__, __FILE__, 1);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2xNCP1(complex<double> *in, complex<double> *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

//  status=projection_DNP1xNCP1_template(in, out,  mesh);
   check_error(-1, "projection not yet implemented", __LINE__, __FILE__, 1);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP0_template(T *in, T *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on LGP0 nodes */
{
  int n,status;

  switch(discretisation) {
    case LGP0:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.ntriangles; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case LGP1:
/**----------------------------------------------------------------------
      */
      status=projection_LGP1xLGP0_template(in, out, mesh);
      break;

    case LGP2:
/**----------------------------------------------------------------------
      */
      status=projection_LGP2xLGP0_template(in, out, mesh);
      break;

    case NCP1:
/**----------------------------------------------------------------------
      */
      status=projection_NCP1xLGP0_template(in, out, mesh);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0(float *in, float *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on LGP0 nodes */
{
  int status;
  status=projection_LGP0_template(in, out, mesh, discretisation);
  return(status);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0(double *in, double *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on LGP0 nodes */
{
  int n,status;

  switch(discretisation) {
    case LGP0:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.ntriangles; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case LGP1:
/**----------------------------------------------------------------------
      */
      status=projection_LGP1xLGP0(in, out, mesh);
      break;

    case LGP2:
/**----------------------------------------------------------------------
      */
      status=projection_LGP2xLGP0(in, out, mesh);
      break;

    case NCP1:
/**----------------------------------------------------------------------
      */
      status=projection_NCP1xLGP0(in, out, mesh);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP0(complex<double> *in, complex<double> *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double complex precision vector projection on LGP0 nodes */
{
  int n,status;

  switch(discretisation) {
    case LGP0:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.ntriangles; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case LGP1:
/**----------------------------------------------------------------------
      */
      status=projection_LGP1xLGP0(in, out, mesh);
      break;

    case LGP2:
/**----------------------------------------------------------------------
      */
//      status=projection_LGP2xLGP0(in, out, mesh);
      check_error(-1, "discretisation not implemented yet", __LINE__, __FILE__, 1);
      break;

    case NCP1:
/**----------------------------------------------------------------------
      */
//      status=projection_NCP1xLGP0(in, out, mesh);
      check_error(-1, "discretisation not implemented yet", __LINE__, __FILE__, 1);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int projection_LGP1_template(T *in, T *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on LGP1 nodes */
{
  int n,status;
  
  if(mesh.LGP1descriptor.massmatrix.packed==0) {
    check_error(-1, "LGP1 descriptor mass matrix not initialised", __LINE__, __FILE__, 1);
    }

  switch(discretisation) {
    case LGP0:
/**----------------------------------------------------------------------
      */
      status=projection_LGP0xLGP1(in, out, mesh);
      break;

    case LGP1:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.nvtxs; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case DGP1:
/**----------------------------------------------------------------------
      */
      status=projection_DGP1xLGP1(in, out, mesh);
      break;

    case LGP2:
/**----------------------------------------------------------------------
      */
      status=projection_LGP2xLGP1(in, out, mesh);
      break;

    case NCP1:
/**----------------------------------------------------------------------
      */
      status=projection_NCP1xLGP1(in, out, mesh);
      break;

    case DNP1:
/**----------------------------------------------------------------------
      */
      status=projection_DNP1xLGP1(in, out, mesh);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int projection_LGP1(double *in, double *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=projection_LGP1_template(in, out, mesh, discretisation);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int projection_LGP1(complex<double> *in, complex<double> *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=projection_LGP1_template(in, out, mesh, discretisation);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP1(float *in, float *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on LGP1 nodes */
{
  int status;
  
  status=projection_LGP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_DGP1_template(T *in, T *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on DGP1 nodes */
{
  int n,status;

  switch(discretisation) {
    case LGP0:
/**----------------------------------------------------------------------
      */
      check_error(-1, "discretisation not implemented", __LINE__, __FILE__, 1);
      status=projection_LGP0xDGP1(in, out, mesh);
      break;

    case LGP1:
/**----------------------------------------------------------------------
      */
      status=projection_LGP1xDGP1(in, out, mesh);
//       check_error(-1, "discretisation not implemented", __LINE__, __FILE__, 1);
//       if(in==out) return(0);
//       for(m = 0; m < mesh.ntriangles; m++) {
//         for(k=0;k<3;k++) {
//           n=mesh.DGP1descriptor.NIbE[m][k];
//           out[n] = in[mesh.triangles[m].vertex[k]];
//           }
//         }
      status=0;
      break;

    case DGP1:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.DGP1descriptor.nnodes; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case NCP1:
/**----------------------------------------------------------------------
      */
      check_error(-1, "discretisation not implemented", __LINE__, __FILE__, 1);
      break;

    case LGP2:
/**----------------------------------------------------------------------
      */
      status=projection_LGP2xDGP1(in, out, mesh);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1(float *in, float *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on DGP1 nodes */
{
  int status;
  
  status=projection_DGP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1(double *in, double *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on DGP1 nodes */
{
  int status;
  
  status=projection_DGP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DGP1(complex<double> *in, complex<double> *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on DGP1 nodes */
{
  int status;
  
  status=projection_DGP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_LGP2_template(T *in, T *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection onto LGP2 nodes */
{
  int n,status;

  switch(discretisation) {
    case LGP0:
      status=projection_LGP0xLGP2(in, out, mesh);
      break;

    case LGP1:
      status=projection_LGP1xLGP2(in, out, mesh);
      break;

    case DGP1:
      status=projection_DGP1xLGP2(in, out, mesh);
      break;

    case NCP1:
      status=projection_NCP1xLGP2(in, out, mesh);
      break;

    case DNP1:
      status=projection_DNP1xLGP2(in, out, mesh);
      break;

    case LGP2:
      if(in==out) return(0);
      for(n = 0; n < mesh.LGP2descriptor.nnodes; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case IPG:
      status=projection_IPGxLGP2(in, out, mesh);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2(float *in, float *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on LGP2 nodes */
{
  int status;
  
  status=projection_LGP2_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2(double *in, double *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on LGP2 nodes */
{
  int status;
  
  status=projection_LGP2_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_LGP2(complex<double> *in, complex<double> *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on LGP2 nodes */
{
  int status;
  
  status=projection_LGP2_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_NCP1_template(T *in, T *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;

  switch(discretisation) {
    case LGP0:
/**----------------------------------------------------------------------
      */
      status=projection_LGP0xNCP1(in, out, mesh);
      break;

    case LGP1:
/**----------------------------------------------------------------------
      */
      status=projection_LGP1xNCP1(in, out, mesh);
      break;

    case DGP1:
/**----------------------------------------------------------------------
      */
      check_error(-1, "discretisation not implemented yet", __LINE__, __FILE__, 1);
      break;

    case NCP1:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.nedges; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case DNP1:
/**----------------------------------------------------------------------
      */
      status=projection_DNP1xNCP1(in, out, mesh);
      break;


    case LGP2:
/**----------------------------------------------------------------------
      */
      status=projection_LGP2xNCP1(in, out, mesh);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1(float *in, float *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_NCP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1(double *in, double *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_NCP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_NCP1(complex<double> *in, complex<double> *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_NCP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_DNP1_template(T *in, T *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on DNP1 nodes */
{
  int n,status;

  switch(discretisation) {
    case LGP0:
/**----------------------------------------------------------------------
      */
      check_error(-1, "discretisation not implemented yet", __LINE__, __FILE__, 1);
      break;

    case LGP1:
/**----------------------------------------------------------------------
      */
      status=projection_LGP1xDNP1(in, out, mesh);
      break;

    case DGP1:
/**----------------------------------------------------------------------
      */
      check_error(-1, "discretisation not implemented yet", __LINE__, __FILE__, 1);
      break;

    case NCP1:
/**----------------------------------------------------------------------
      */
      status=projection_NCP1xDNP1(in, out, mesh);
      break;

    case DNP1:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.DNP1descriptor.nnodes; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case LGP2:
/**----------------------------------------------------------------------
      */
      status=projection_LGP2xDNP1(in, out, mesh);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1(float *in, float *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_DNP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1(double *in, double *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_DNP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_DNP1(complex<double> *in, complex<double> *out, mesh_t & mesh, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=projection_DNP1_template(in, out, mesh, discretisation);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_IPG_template(T *in, T *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status=0;
  T *p;
  gauss_t gauss;
  
  discretisation_t src_descriptor=get_descriptor(mesh,source);
  discretisation_t dst_descriptor=mesh.IPGdescriptor;

  gauss.init(mesh.IPGdescriptor.type);
  
  for(n = 0; n < dst_descriptor.nnodes; n++)
    out[n] = 0;

  if(in==0) return(-1);
  
  p=new T[src_descriptor.nnpe];
  
  for (m=0;m<mesh.ntriangles;m++) {
    for(k=0;k<src_descriptor.nnpe;k++) {
      n=src_descriptor.NIbE[m][k];
      p[k]=in[n];
      }
    for(k=0;k<dst_descriptor.nnpe;k++) {
      n=dst_descriptor.NIbE[m][k];
      out[n] = gauss.interpolate(p,k,src_descriptor.type);
      }
    }

  gauss.destroy();
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_IPG(float *in, float *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_IPG_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_IPG(double *in, double *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_IPG_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_IPG(complex<double> *in, complex<double> *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_IPG_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_CQP1xCQP0_template(T *in, T *out, mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n, status;
  
  for(n = 0; n < mesh.CQP0descriptor.nnodes; n++) {
    out[n] = 0;
    for(k=0;k<4;k++) {
      m=mesh.CQP0descriptor.nodes[n].vtces[k];
      out[n]+=in[m];
      }
    out[n] /=4.0;;
    }

//   for (m=0;m<mesh.nquadrangles;m++) {
//     for(k=0;k<4;k++) {
//       n=mesh.quadrangles[m].vertex[k];
//       p[k]=in[n];
//       }
//     sum=0;
//     for(k=0;k<3;k++) {
//       sum+=p[k];
//       }
//     out[m]=sum/3.0;
//     }
     
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_CQP0_template(T *in, T *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  switch(source) {

    case CQP0:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.CQP0descriptor.nnodes; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case LGP1: /// HERE !!!
    case CQP1:
/**----------------------------------------------------------------------
      */
      status=projection_CQP1xCQP0_template(in, out, mesh);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQP0(float *in, float *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQP0_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQP0(double *in, double *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQP0_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQP0(complex<double> *in, complex<double> *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQP0_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_CQP1_template(T *in, T *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  T* count;
  
  switch(source) {

    case CQP0:
/**----------------------------------------------------------------------
      */
      status=0; /// HERE !!!
      break;

    case LGP1: /// HERE !!!
    case CQP1:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.CQP1descriptor.nnodes; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    case CQN1:
/**----------------------------------------------------------------------
      */
      check_error(-1, "CQN1xCQP1 projection to be implemented properly", __LINE__, __FILE__, 0);
      count=new T[mesh.CQP1descriptor.nnodes];
      for(n = 0; n < mesh.CQP1descriptor.nnodes; n++) count[n]=0;
      for(n = 0; n < mesh.CQP1descriptor.nnodes; n++) out[n]=0;
      for(n = 0; n < mesh.CQN1descriptor.nnodes; n++) {
        size_t n1=mesh.edges[n].extremity[0];
        size_t n2=mesh.edges[n].extremity[1];
        out[n1] += in[n];
        count[n1]+=1;
        out[n2] += in[n];
        count[n2]+=1;
        }
      for(n = 0; n < mesh.CQP1descriptor.nnodes; n++) out[n]/=count[n];
      status=0;
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQP1(float *in, float *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQP1_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQP1(double *in, double *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQP1_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQP1(complex<double> *in, complex<double> *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQP1_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int projection_CQN1_template(T *in, T *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  switch(source) {

    case CQP0:
/**----------------------------------------------------------------------
      */
      status=0; /// HERE !!!
      break;

    case LGP1: /// HERE !!!
    case CQP1:
/**----------------------------------------------------------------------
      */
      for(n = 0; n < mesh.CQN1descriptor.nnodes; n++) {
        size_t n1=mesh.edges[n].extremity[0];
        size_t n2=mesh.edges[n].extremity[1];
        out[n] = (T) 0.5*(in[n1]+in[n2]);
        }
      status=0;
      break;

    case CQN1:
/**----------------------------------------------------------------------
      */
      if(in==out) return(0);
      for(n = 0; n < mesh.CQN1descriptor.nnodes; n++) {
        out[n] = in[n];
        }
      status=0;
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQN1(float *in, float *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQN1_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQN1(double *in, double *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQN1_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQN1(complex<float> *in, complex<float> *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQN1_template(in, out, mesh,source));
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQN1(complex<double> *in, complex<double> *out, mesh_t & mesh, int source)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    return(projection_CQN1_template(in, out, mesh,source));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int fe_projection_template(mesh_t & mesh, T *in, int source, T *out, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------
  Perform double precision vector projection on target nodes */
{
  int status;

  switch(target) {
    case LGP0:
/**----------------------------------------------------------------------
      */
      status=projection_LGP0_template(in, out, mesh, source);
      break;

    case LGP1:
/**----------------------------------------------------------------------
      */
      status=projection_LGP1_template(in, out, mesh, source);
      break;

    case DGP1:
/**----------------------------------------------------------------------
      */
      status=projection_DGP1_template(in, out, mesh, source);
      break;

    case LGP2:
/**----------------------------------------------------------------------
      */
      status=projection_LGP2_template(in, out, mesh, source);
      break;

    case NCP1:
/**----------------------------------------------------------------------
      */
      status=projection_NCP1_template(in, out, mesh, source);
      break;

    case DNP1:
/**----------------------------------------------------------------------
      */
      status=projection_DNP1_template(in, out, mesh, source);
      break;

    case IPG:
/**----------------------------------------------------------------------
      */
      status=projection_IPG_template(in, out, mesh, source);
      break;

    case CQP0:
/**----------------------------------------------------------------------
      */
      status=projection_CQP0_template(in, out, mesh, source);
      break;

    case CQP1:
/**----------------------------------------------------------------------
      */
      status=projection_CQP1_template(in, out, mesh, source);
      break;

    case CQN1:
/**----------------------------------------------------------------------
      */
      status=projection_CQN1_template(in, out, mesh, source);
      break;

    default:
      check_error(-1, "unknown discretisation", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_projection(mesh_t & mesh, float *in, int source, float *out, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=fe_projection_template( mesh, in, source, out,  target);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_projection(mesh_t & mesh, double *in, int source, double *out, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=fe_projection_template( mesh, in, source, out,  target);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_projection(mesh_t & mesh, complex<float> *in, int source, complex<float> *out, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=fe_projection_template( mesh, in, source, out,  target);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_projection(mesh_t & mesh, complex<double> *in, int source, complex<double> *out, int target)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=fe_projection_template( mesh, in, source, out,  target);
  
  return status;
}


/**----------------------------------------------------------------------------
  alternative u compuation from U; disabled */
//   tmp=new double[u_descriptor.nnodes];
//   for(n = 0; n < udim; n++) {
//     tmp[n]=0;
//     }
//
//   for (m=0;m<mesh.ntriangles;m++) {
//     double pseudo=mesh.triangles[m].Area;
//     for(k=0;k<u_descriptor.nnpe;k++) {
//       n=u_descriptor.NIbE[m][k];
//       p[k]=state[next].U[n];
//       }
//     for(k=0;k<z_descriptor.nnpe;k++) {
//       n=z_descriptor.NIbE[m][k];
//       q[k]=1./state[next].H[n];
//       }
//     double C=mesh.triangles[m].cosinus;
//     for(k=0;k<u_descriptor.nnpe;k++) {
//       n=u_descriptor.NIbE[m][k];
//       tmp[n] += 2.0 * pseudo * C * fe_integraleLGP2xNCP1xNCP1_2D(q,p,k);
//       }
//     }
//   status = LinearSystem_solve(u_descriptor.massmatrix,tmp);
//   for(n = 0; n < udim; n++) {
//     state[next].u[n]=tmp[n];
//     }
//
//   for(n = 0; n < udim; n++) {
//     tmp[n]=0;
//     }
//
//   for (m=0;m<mesh.ntriangles;m++) {
//     double pseudo=mesh.triangles[m].Area;
//     for(k=0;k<u_descriptor.nnpe;k++) {
//       n=u_descriptor.NIbE[m][k];
//       p[k]=state[next].V[n];
//       }
//     for(k=0;k<z_descriptor.nnpe;k++) {
//       n=z_descriptor.NIbE[m][k];
//       q[k]=1./state[next].H[n];
//       }
//     double C=mesh.triangles[m].cosinus;
//     for(k=0;k<u_descriptor.nnpe;k++) {
//       n=u_descriptor.NIbE[m][k];
//       tmp[n] += 2.0 * pseudo * C * fe_integraleLGP2xNCP1xNCP1_2D(q,p,k);
//       }
//     }
//   status = LinearSystem_solve(u_descriptor.massmatrix,tmp);
//   for(n = 0; n < udim; n++) {
//     state[next].v[n]=tmp[n];
//     }
//
//   delete[] tmp;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int projection_CQP1xCQP0(mesh_t & mesh, double *in, double *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,m,n,status=0;
//   double pseudo,C;
  discretisation_t descriptor=mesh.CQP0descriptor;

  if(in==0) return(-1);

  for(n = 0; n < descriptor.nnodes; n++)
    out[n] = 0;

  for (m=0;m<mesh.nquadrangles;m++) {
//    pseudo=mesh.quadrangles[m].Area;
    for(k=0;k<4;k++) {
      n=mesh.quadrangles[m].vertex[k];
      out[m]+=in[n];
      }
    out[m]/=4.0;
//    C=mesh.quadrangles[m].cosinus;
    }

  return(status);
}

