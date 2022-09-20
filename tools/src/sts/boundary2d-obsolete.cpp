
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

#define LAPACKF_ 1

#include "tugo-prototypes.h"

// void Boundary_tide_obsolete(double, z_spec *, int, float *, float *, float *);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Hindex(int row, int col, int hbw)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;

  k = col * (3 * hbw + 1) + row - col + 2 * hbw;
  if(row - col + 2 * hbw < hbw) {
    printf("warning\n");
    }
  if(row - col + 2 * hbw > 3 * hbw + 1) {
    printf("warning\n");
    }
  return (k);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int Vindex(int i, int j, int hbw)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;
  k = i * hbw + (i - j) + hbw + 1;
  return (k);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dBoundaryMatrix(mesh_t mesh, int level, double *A, int *pivot, int hbw)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------------

  Variational formulation of Dirichlet boundary conditions

------------------------------------------------------------------------------*/
{
  int k, m, n, n1, n2;
  int bw, bm, neq, status = -1;
  int i, j;
  triangle_t *elt;
  int nndes, nelts;
  double *in, *out;

  bw = 2 * hbw + 1;
  bm = 3 * hbw + 1;

  for(k = 0; k < gOpenBCs[level].nnodes * bm; k++)
    A[k] = 0;

/*----------------------------------------------------------------------
  */

  for(k = 0; k < gOpenBCs[level].nedges; k++) {
    n  = gOpenBCs[level].edges[k];
    n1 = mesh.edges[n].extremity[0];
    n2 = mesh.edges[n].extremity[1];
    }

/*----------------------------------------------------------------------
  build linear element mass matrix*/

  for(k = 0; k < gOpenBCs[level].nnodes; k++) {
    m = gOpenBCs[level].r_edge[k];
    n = gOpenBCs[level].r_node[k];
    if(n != -1) {
      A[Hindex(k, n, hbw)] += 1. / 6. * mesh.edges[m].L;
      A[Hindex(k, k, hbw)] += 1. / 3. * mesh.edges[m].L;
      }
    m = gOpenBCs[level].l_edge[k];
    n = gOpenBCs[level].l_node[k];
    if(n != -1) {
      A[Hindex(k, n, hbw)] += 1. / 6. * mesh.edges[m].L;
      A[Hindex(k, k, hbw)] += 1. / 3. * mesh.edges[m].L;
      }
    }
//   for(k = 0; k < gOpenBCs[level].nnodes; k++) {
//     A[Hindex(k, k, hbw)] = 1;
//     }

  neq = gOpenBCs[level].nnodes;

#ifdef LINPACKF_
  dgbfa_(A, &bm, &neq, &hbw, &hbw, pivot, &status);
#elif LINPACKF
  dgbfa(A, &bm, &neq, &hbw, &hbw, pivot, &status);
#elif LINPACKC
  dgbfa(A, bm, neq, hbw, hbw, pivot, &status);
#elif  LAPACKF_
  dgbtrf_(&neq, &neq, &hbw, &hbw, A, &bm, pivot, &status);
#elif  LAPACKC
  dgbtrf(neq, neq, hbw, hbw, A, bm, pivot, &status);
#elif ATLAS
// "102" changed to "CblasColMajor"
//  status = clapack_dgbtrf(CblasColMajor, neq, neq, A, neq, pivot);
#endif

  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void dBoundarySolver(mesh_t mesh, int level, double *A, int *pivot, int hbw, double *in, double *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**-----------------------------------------------------------------------------

  Variational formulation of Dirichlet boundary conditions

------------------------------------------------------------------------------*/
  int i, j, k, l, m, n, n1, n2, status;
  edge_t *edges;

  int nelts, nedges, nndes;
  double p[3], q[3];
  triangle_t elt;
  int neq, bandmax, bandwidth, ml, mu, job,nrhs=1;

  ml = hbw;
  mu = hbw;

  bandmax = 2 * ml + mu + 1;
  bandwidth = ml + mu + 1;

  for(n = 0; n < gOpenBCs[level].nnodes; n++)
    out[n] = 0;

  for(k = 0; k < gOpenBCs[level].nnodes; k++) {
    m = gOpenBCs[level].r_edge[k];
    n = gOpenBCs[level].r_node[k];
    if(n != -1) {
      out[k] += 1. / 6. * mesh.edges[m].L * in[n];
      out[k] += 1. / 3. * mesh.edges[m].L * in[k];
      }
    m = gOpenBCs[level].l_edge[k];
    n = gOpenBCs[level].l_node[k];
    if(n != -1) {
      out[k] += 1. / 6. * mesh.edges[m].L * in[n];
      out[k] += 1. / 3. * mesh.edges[m].L * in[k];
      }
    }

  neq = gOpenBCs[level].nnodes;

  job = 0;

#ifdef LINPACKF_
  dgbsl_(A, &bandmax, &neq, &hbw, &hbw, pivot, out, &job);
#elif LINPACKF
  dgbsl(A, &bandmax, &neq, &hbw, &hbw, pivot, out, &job);
#elif LINPACKC
  dgbsl(A, bandmax, neq, hbw, hbw, pivot, out, job);
#elif  LAPACKF_
  char task='N';
  dgbtrs_(&task, &neq, &hbw, &hbw, &nrhs, A, &bandmax, pivot, out, &neq, &status);
  if(status!=0) {
    check_error(-1, "dgbtrs (double band matrix) solver failed", __LINE__, __FILE__, 1);
    }
#elif  LAPACKC
  char task='N';
  dgbtrs_(task, neq, hbw, hbw, nrhs, A, bandmax, pivot, out, neq, &status);
  if(status!=0) {
    check_error(-1, "dgbtrs (double band matrix) solver failed", __LINE__, __FILE__, 1);
    }
#elif ATLAS
// "102" changed to "CblasColMajor"
//  status = clapack_dgetrs(CblasColMajor, CblasNoTrans, neq, nrhs, A, neq, pivot, out, neq);
#else
  check_error(-1, "no linear algebra library available", __LINE__, __FILE__, 1);
#endif

}



