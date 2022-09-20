/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2017

  Unstructured Ocean Grid initiative

Developers:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
Contributors:
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France
  Frédéric Dupont    Canada environement, Canada

***************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <time.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "solvers-functions.h"
#include "linear-gmres.h"

#if 0

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class C> int operation(HYPERMATRIX_t<C> & M, HYPERMATRIX_t<C> & P, C *X, C *Y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t size=M.ordering->nrows;
  
  C *r=new C[size];
  
  status=matrix_operation(M, X, r, 1);
  status=matrix_operation(P, r, Y, 1);
  
  delete[] r;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class C> void arnoldi(HYPERMATRIX_t<C> & A, HYPERMATRIX_t<C> & P, C *Q, vector_template_t<C> & q, C *h, int k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Arnoldi function: 
  
  add the next (k+1 th) orthogonal vector to the Krylov vector basis (Q), i.e. q
  
  add the next (k th) column of the Hessenberg matrix (H), i.e. h
  
  strong MPI constraints
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  
  size_t nrows=A.neq;
  
//   status=matrix_operation(A, &Q[k*nrows], q.x, 1);
  status=operation(A, P, &Q[k*nrows], q.x);
  
/*------------------------------------------------------------------------------
  WARNING : MPI sync of q will be needed here (partitions boundary nodes) */  
  
  for (int i=0; i<k+1; i++) {
    vector_template_t<C> p;
/*------------------------------------------------------------------------------
    p is i th columns of Q, i.e. i th vetor of Krylov basis */
    p.duplicate(&Q[i*nrows], nrows);
    h[i]=q*p;
/*------------------------------------------------------------------------------
    WARNING : MPI sum will be needed to achieve global h[i] */  
    p*=h[i];
/*------------------------------------------------------------------------------
    WARNING : MPI sync of q will be needed here (partitions boundary nodes) */  
    q-=p;
    }

  h[k+1] = norm(q);
  q/=h[k+1];
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void givens_rotation(double v1, double v2, double & cs, double & sn)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Calculate the Given rotation matrix
  
  no MPI constraint
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  if (v1==0) {
    cs = 0;
    sn = 1;
    }
  else {
    double t=sqrt(v1*v1+v2*v2);
    cs = abs(v1) / t;
    sn = cs * v2 / v1;
//     cs = v1 / t;
//     sn = v2 / t;
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void givens_rotation(complex<double> v1, complex<double> v2, complex<double> & cs, complex<double> & sn)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Calculate the Given rotation matrix
  
  no MPI constraint
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double tmp=abs(v1);
  
  if (tmp==0) {
    cs = 0;
    sn = 1;
    }
  else {
    complex<double> t=sqrt(tmp*tmp+v2*v2);
    cs = abs(v1) / t;
    sn = cs * v2 / v1;
//     cs = v1 / t;
//     sn = v2 / t;
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void apply_givens_rotation(double *h, double *cs, double *sn, int k, double & cs_k,double & sn_k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Applying Givens Rotation to H col
  
  no MPI constraint
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  
/*------------------------------------------------------------------------------
  apply for i-th column */
  for (int i=0; i<k; i++) {
    double temp;
    temp   =  cs[i]*h[i] + sn[i]*h[i+1];
    h[i+1] = -sn[i]*h[i] + cs[i]*h[i+1];
    h[i]   = temp;
    }
  
/*------------------------------------------------------------------------------
  update the next sin cos values for rotation */
  givens_rotation( h[k], h[k+1], cs_k, sn_k);
  
/*------------------------------------------------------------------------------
  eliminate H(i+1,i) */
  h[k]   = cs_k*h[k] + sn_k*h[k+1];
  h[k+1] = 0.0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void apply_givens_rotation(complex<double> *h, complex<double> *cs, complex<double> *sn, int k, complex<double> & cs_k,complex<double> & sn_k)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Applying Givens Rotation to H col
  
  no MPI constraint
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  
/*------------------------------------------------------------------------------
  apply for i-th column */
  for (int i=0; i<k; i++) {
    complex<double> temp;
    temp   =  conj(cs[i]*h[i]) + conj(sn[i])*h[i+1];
    h[i+1] = -sn[i]*h[i] + cs[i]*h[i+1];
    h[i]   = temp;
    }
  
/*------------------------------------------------------------------------------
  update the next sin cos values for rotation */
  givens_rotation( h[k], h[k+1], cs_k, sn_k);
  
/*------------------------------------------------------------------------------
  eliminate H(i+1,i) */
  h[k]   = conj(cs_k)*h[k] + conj(sn_k)*h[k+1];
  h[k+1] = 0.0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class C> int UpperTriangular_solver(C *U, int size, C *B, C *x, int n,int & singular)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  solver:
    solve the linear equation Ux = B for x, where U is an upper
    triangular matrix.                                      
 
  size is the allocated row dimension (potentially greater than n)
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
   int i, k;

   singular=-1;
   
   for (k = n-1; k >= 0; k--) {
     size_t pos=k*size+k;
     C c=U[pos];
     if (c == 0.0) {
       singular=k;
       return -1;           // The matrix U is singular
       }
     x[k] = B[k];
     for (i = k + 1; i < n; i++) {
       size_t pos=i*size+k;
       x[k] -= x[i] * U[pos];
       }
     x[k] /= c;
     }

   return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <class C> int gmres_solution(gmres_t<C> & gmres, C *x, C *H, C *Q, vector_template_t<C>  & g, int niter, int m, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  solve for the least square minimisation 
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;

  C *y=gmres.y;
  int singular;
  
  status=UpperTriangular_solver(H, m+1, g.x, y, niter, singular);
  if(status!=0) {
    printf("UpperTriangular_solver failed at line %d\n",singular);
    return(-1);
    }
    
/*------------------------------------------------------------------------------
  update x vector */
  for(int i=0; i<n;i++) {
    C sum=0.0;
    for(int k=0; k<niter; k++) {
      sum+=Q[k*n+i]*y[k];
      }
    x[i]+=sum;
    }  

#if 0
  for(int k=0; k<niter; k++) {
    dvector_t q;
    q.alias(&Q[n*k], n);
    printf("%s : k %3d y=%g q=%g\n", __func__, k, y[k], norm(q));
    }
#endif

  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//   int gmres_solver(mesh_t & mesh, hypermatrix_t & A, hypermatrix_t & P, double *rhs, double *x, double *e, int max_iterations, double threshold, gmres_t & gmres) // GMRES IMPLEMENTATION, TEMPORARILY COMMENTED
  template <class C> int gmres_solver_template(HYPERMATRIX_t<C> & A, HYPERMATRIX_t<C> & P, C *rhs, C *x, double *e, int max_iterations, double threshold, gmres_t<C> & gmres)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  GMRES solver
  
  locally low MPI constraints

  inspired from : 
  
    https://en.wikipedia.org/wiki/Generalized_minimal_residual_method
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  size_t n = A.neq;
  int m=max_iterations, kmax=max_iterations-2;
  bool debug=false;
//   static gmres_t gmres;
  int verbose=0;
  
  if(not gmres.initialised) gmres.init(max_iterations, n); 
//   else status=gmres_solution(gmres, x, gmres.H, gmres.Q, gmres.g, gmres.niter, gmres.nitermax, gmres.nunknowns);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

   compute pre-conditioned residual: r=P(b-A*x0)
   
   use entering x as the initial solution approximation vector x0
        
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  vector_template_t<C> &r=gmres.r;
  
  vector_template_t<C> b;
  b.alias(rhs,n);
  
  if(debug) {
    status=vector_checkInf(b.x, n);
    status=vector_checkNaN(b.x, n);
    }
     
  status=operation(A, P, x, r.x);
/*------------------------------------------------------------------------------
  WARNING : MPI sync of r will be needed here (partitions boundary nodes) */  
  
  if(debug) {
    status=vector_checkNaN(x,   n);
    status=vector_checkNaN(r.x, n);
    }
  
  r*=-1.0;
  status=matrix_operation(P, b.x, r.x, 0);
/*------------------------------------------------------------------------------
  WARNING : MPI sync of r will be needed here (partitions boundary nodes) */  

  double b_norm = norm(b);
  
  if(isinf(b_norm)) {
    TRAP_ERR_EXIT(-1, "nan encountered\n");
    }

/*------------------------------------------------------------------------------
  WARNING : MPI sum will be needed for b_norm */

  double r_norm=norm(r);
/*------------------------------------------------------------------------------
  WARNING : MPI sum will be needed for r_norm */

  double error  = r_norm/b_norm;
/*------------------------------------------------------------------------------
  WARNING : MPI sum will be needed for error */  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 initialize the 1D vectors
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  vector_template_t<C> &cs=gmres.cs; // HERE !!! to be fixed
  vector_template_t<C> &sn=gmres.sn; // HERE !!! to be fixed
  vector_template_t<C>  &g =gmres.g;
  dvector_t &e1=gmres.e1;
  
  e1.x[0] = 1.0;
  e[0]=error;
  if(verbose==1) printf("%s : iteration %3d error=%g\n", __func__, 0, error);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  initialize Q and H matrix
 
  actually H is not needed, but R
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  C *Q=gmres.Q;
  C *H=gmres.H;
 
/*------------------------------------------------------------------------------
  initialise first q with normalized residual */
  for(int k=0; k<n; k++) Q[k]=r.x[k]/r_norm;
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  try to use previous Krylov space to improve first guess
  
  buggy for first pass, needs to be reworked and assessed
  
  also, static gmres_t not convenient for subcycling
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
//   status=gmres_solution(gmres, x, H, Q, g, 1, m, n); 
//   status=force_periodic(mesh, gPeriodic[0], x);
//   
//   status=operation(A, P, x, r.x);
//   r*=-1.0;
//   status=matrix_operation(P, b.x, r.x, 0);
//   r_norm=norm(r);
//   error  = r_norm/b_norm;
//   if(verbose==1) printf("%s : iteration %3d error=%g\n", __func__, 0, error);
//   for(int k=0; k<n; k++) Q[k]=r.x[k]/r_norm;

/*------------------------------------------------------------------------------
  WARNING : MPI sync of Q will be needed here (partitions boundary nodes) */  

/*------------------------------------------------------------------------------
  initialise first q with residual norm */
  g.x[0] = r_norm;
  
  bool stop=false;
  for(int k=0; k<m-1; k++) {
    vector_template_t<C> q;
    q.alias(&Q[n*(k+1)], n);
    
/*------------------------------------------------------------------------------
    Arnoldi iteration : add k+1 th orthogonal vector to Krylov space basis */
    arnoldi(A, P, Q, q, &H[k*(m+1)], k);

/*------------------------------------------------------------------------------
    eliminate the last element in H ith row and update the rotation matrix */
    apply_givens_rotation(&H[k*(m+1)], cs.x, sn.x, k, cs.x[k], sn.x[k]); // HERE !!! to be fixed

/*------------------------------------------------------------------------------
    update the residual vector */
    g.x[k+1] = -sn.x[k]*g.x[k]; // HERE !!! to be fixed
    g.x[k]   =  cs.x[k]*g.x[k]; // HERE !!! to be fixed
    
    error=abs(g.x[k+1]) / b_norm;
    
//     printf("%s : iteration %3d error=%g\n", __func__, k, error);
    
/*------------------------------------------------------------------------------
    save the error */
    e[k+1]=error;
    
    if(stop) {
      kmax++;
      break;
      }
//     if(error <= threshold/b_norm) {
    if(error <= 1.e-10) {
      kmax=k;
      stop=true;
      }
    }

  gmres.niter=kmax+1;
  if(error>1.e-10) printf("%s : iteration %3d error=%g\n", __func__, 0, error);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  calculate the solution
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=gmres_solution(gmres, x, H, Q, g, kmax+1, m, n);

//   status=force_periodic(mesh, gPeriodic[0], x); // GMRES IMPLEMENTATION, TEMPORARILY COMMENTED
    
/*------------------------------------------------------------------------------
  WARNING : MPI sync of x will be needed here (partitions boundary nodes) */  
  
//   printf("%s : iteration %3d error=%g\n", __func__, kmax, error);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  verification
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(debug) {
    status=operation(A, P, x, r.x);
    r*=-1.0;
    status=matrix_operation(P, b.x, r.x, 0);
//     status=matrix_operation(A, x, r.x, 1);
//     r*=-1.0;
//     r+=b;
    r_norm = norm(r);
    }
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//   int gmres_solver(mesh_t & mesh, hypermatrix_t & A, hypermatrix_t & P, double *rhs, double *x, double *e, int max_iterations, double threshold, gmres_t & gmres) // GMRES IMPLEMENTATION, TEMPORARILY COMMENTED
  int gmres_solver(hypermatrix_t & A, hypermatrix_t & P, double *rhs, double *x, double *e, int max_iterations, double threshold, gmres_t<double> & gmres)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  wrapper
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;

  status=gmres_solver_template(A, P, rhs, x, e, max_iterations, threshold, gmres);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//   int gmres_solver(mesh_t & mesh, hypermatrix_t & A, hypermatrix_t & P, double *rhs, double *x, double *e, int max_iterations, double threshold, gmres_t & gmres) // GMRES IMPLEMENTATION, TEMPORARILY COMMENTED
  int gmres_solver(hyperzmatrix_t & A, hyperzmatrix_t & P, complex<double> *rhs, complex<double> *x, double *e, int max_iterations, double threshold, gmres_t< complex<double> > & gmres)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  wrapper
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;

  status=gmres_solver_template(A, P, rhs, x, e, max_iterations, threshold, gmres);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmres(hypermatrix_t& M, complex<double>  *x, complex<double>  *y, int init, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  wrapper
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  
  status=matrix_operation(M, x, y, init);
  
/*------------------------------------------------------------------------------
  complete solution by node exchanges for unsolved (overlapped) unknowns*/
  if(M.mpi_exchange==0) TRAP_ERR_EXIT(-1, "%s : mpi_exchange not initialised (matrix %s)\n",__func__,M.name.c_str());

//   status=exchange_atom(*M.mpi_exchange, gCPU_ID, gnCPUs, y, y, debug);
  
  
  
  return(status);
  
}

#endif
