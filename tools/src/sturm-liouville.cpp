
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2011

  Unstructured Ocean Grid initiative

Contributors:
 
  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  Clement Mayet      LEGOS, Toulouse, France
  Yves Soufflet      LEGOS, Toulouse, France
  Damien Allain      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"
#include "fe.h"
#include "map.h"
#include "grd.h"
#include "matrix.h"
#include "netcdf-proto.h"

#include "solverlib.h"

extern "C" double ddot_  (int *,double *,int *,double *,int *);

extern "C" void dgeev_(char *, char *, int *, double *, int *, double *, double *, double *, int *, double *, int *, double *, int *, int *);
extern "C" void dggev_(char *, char *, int *, double *, int *, double *, int *, double *, double *, double *, double *, int *, double *, int *, double *, int *, int *);

extern double water_density(double t, double s, double z, int option=2);
extern double water_density_check(int option);
extern double physic_theta(double s4, double t04, double p04, double pr);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int check_Wmodes(int nlevels, double *z, double **modes, int nmodes, double *rho, double *N, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k;
  int nlayers=nlevels-1;
  double norm,scalar_product;
  
  double *dz=new double[nlayers];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int k=0;k<nlayers;k++) {
    dz[k] = z[k+1]-z[k];
    }
  
/**----------------------------------------------------------------------------
  normalize eignevectors to unity for euclidian norm*/
  for(i=0;i<nmodes;i++) {
    norm=0;
    for (k=0;k<nlevels;k++) {
      norm+=modes[i][k]*modes[i][k];
      }
    if(norm==0) {
      if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
      continue;
      }
    norm=sqrt(norm);
    if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
    for (k=0;k<nlevels;k++) {
      modes[i][k]/=norm;
      }
    }
    
/**----------------------------------------------------------------------------
  */
//   for(k=0;k<nlayers;k++) {
//     printf("mode 0, f'= %lf\n",modes[0][k+1]-modes[0][k]);
//     }
//   for(k=0;k<nlayers-1;k++) {
//     printf("mode 0, f''= %12.9lf\n",modes[0][k+2]+modes[0][k]-2*modes[0][k+1]);
//     }
    
/**----------------------------------------------------------------------------
  check orthonormality of eigenvetors basis for euclidian norm*/
  for(i=0;i<nmodes;i++) {
    for(j=i;j<nmodes;j++) {
      scalar_product=0;
      for (k=0;k<nlevels;k++) {
        scalar_product+=modes[i][k]*modes[j][k];
        }
      if(verbose==1) printf("modes %d %d, norm =%lf\n",i,j,scalar_product);
      }
    }
   
/**----------------------------------------------------------------------------
  normalize eignevectors to unity for integrale norm*/
  for(i=0;i<nmodes;i++) {
    norm=0;
    for (k=0;k<nlayers;k++) {
      norm+=0.25*dz[k]*(modes[i][k]+modes[i][k+1])*(modes[i][k]+modes[i][k+1]);
      }
    if(norm==0) {
      if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
      continue;
      }
    norm=sqrt(norm);
    if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
    for (k=0;k<nlevels;k++) {
      modes[i][k]/=norm;
      }
    }
  
/**----------------------------------------------------------------------------
  check orthonormality of eigenvectors basis for integrale norm*/
  for(i=0;i<nmodes;i++) {
    for(j=i;j<nmodes;j++) {
      scalar_product=0;
      for (k=0;k<nlevels-1;k++) {
        scalar_product+=0.25*dz[k]*(modes[i][k]+modes[i][k+1])*(modes[j][k]+modes[j][k+1]);
        }
      if(verbose==1) printf("modes %d %d, norm =%lf\n",i,j,scalar_product);
      }
    }
    
  for(i=0;i<nmodes;i++) {
    for(j=i;j<nmodes;j++) {
      scalar_product=0;
      for (k=0;k<nlevels-1;k++) {
        scalar_product+=0.5*dz[k]*(modes[i][k]*modes[j][k]+modes[i][k+1]*modes[j][k+1]);
        }
      if(verbose==1) printf("#modes %d %d, norm =%lf\n",i,j,scalar_product);
      }
    }
    
//   verbose=1;
  
/**----------------------------------------------------------------------------
  normalize eignevectors to unity for integrale norm with weight function*/
//   for(i=0;i<nmodes;i++) {
//     norm=0;
//     double p[2], q[2], r[2], s[2];
//     for (k=0;k<nlevels-1;k++) {
//       p[0]=N[k]*N[k];
//       p[1]=N[k+1]*N[k+1];
//       norm+=rho[k]*dz[k]*fe_integraleLGP1xLGP1xLGP1_1D(p, &modes[i][k], &modes[i][k]);
//       }
//     if(norm==0) {
//       if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
//       continue;
//       }
//     norm=sqrt(norm);
//     if(verbose==1) printf("W-mode %d, norm =%lf\n",i,norm);
//     for (k=0;k<nlevels;k++) {
//       modes[i][k]/=norm;
//       }
//     }
//     
//   for(i=0;i<nmodes;i++) {
//     for(j=i;j<nmodes;j++) {
//       scalar_product=0;
//       double p[2], q[2], r[2], s[2];
//       for (k=0;k<nlevels-1;k++) {
//         p[0]=N[k]*N[k];
//         p[1]=N[k+1]*N[k+1];
//         scalar_product+=rho[k]*dz[k]*fe_integraleLGP1xLGP1xLGP1_1D(p, &modes[i][k], &modes[j][k]);
//         }
//       if(verbose==1) printf("W-modes %d %d, scalar product =%lf\n",i,j,scalar_product);
//       }
//     }
    
  for(i=0;i<nmodes;i++) {
    norm=0;
    double p[2];
    for (k=1;k<nlevels-1;k++) {
      p[0]=N[k]*N[k];
      p[1]=N[k+1]*N[k+1];
      norm+=rho[k]*N[k]*N[k]*(dz[k-1]+dz[k])*modes[i][k]*modes[i][k]/2.0;
      }
    k=0;
    norm+=rho[k]*N[k]*N[k]*dz[k]*modes[i][k]*modes[i][k]/2.0;
    k=nlevels-1;
    norm+=rho[k-1]*N[k]*N[k]*dz[k-1]*modes[i][k]*modes[i][k]/2.0;
    if(norm==0) {
      if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
      continue;
      }
    norm=sqrt(norm);
    if(verbose==1) printf("W-mode %d, norm =%lf\n",i,norm);
    for (k=0;k<nlevels;k++) {
      modes[i][k]/=norm;
      }
    }
    
  for(i=0;i<nmodes;i++) {
    for(j=i;j<nmodes;j++) {
      scalar_product=0;
      double p[2];
      for (k=1;k<nlevels-1;k++) {
        p[0]=N[k]*N[k];
        p[1]=N[k+1]*N[k+1];
        scalar_product+=rho[k]*N[k]*N[k]*(dz[k-1]+dz[k])*modes[i][k]*modes[j][k]/2.0;
        }
      k=0;
      scalar_product+=rho[k]*N[k]*N[k]*dz[k]*modes[i][k]*modes[j][k]/2.0;
      k=nlevels-1;
      scalar_product+=rho[k-1]*N[k]*N[k]*dz[k-1]*modes[i][k]*modes[j][k]/2.0;
      if(verbose==1) printf("W-modes %d %d, scalar product =%lf\n",i,j,scalar_product);
      }
    }
   
  return 0;
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int normalize_Zmodes(int nlevels, double *z, double **modes, int nmodes, double *rho, double *N, int verbose)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int i, j, status;
//   int nlayers=nlevels-1;
//   double norm,scalar_product;
//   
//   double *dz=new double[nlayers];
// 
// /**----------------------------------------------------------------------------
//   z is the level depth; dz is layer thickness*/
//   for (int k=0;k<nlayers;k++) {
//     dz[k] = z[k+1]-z[k];
//     }
//   
// //   verbose=1;
// 
// /**----------------------------------------------------------------------------
//   normalize eignevectors to unity for integrale norm with weight function*/
//   for(i=0;i<nmodes;i++) {
//     norm=0;
//     double p[2], q[2], r[2], s[2];
//     for (int k=0;k<nlevels-1;k++) {
//       p[0]=N[k]*N[k];
//       p[1]=N[k+1]*N[k+1];
//       norm+=rho[k]*dz[k]*fe_integraleLGP1xLGP1xLGP1_1D(p, &modes[i][k], &modes[i][k]);
//       }
//     norm=sqrt(norm);
//     if(verbose==1) printf("W-mode %d, norm =%lf\n",i,norm);
//     for (int k=0;k<nlevels;k++) {
//       modes[i][k]/=norm;
//       }
//     }
//     
//   for(i=0;i<nmodes;i++) {
//     for(j=i;j<nmodes;j++) {
//       scalar_product=0;
//       double p[2], q[2], r[2], s[2];
//       for (int k=0;k<nlevels-1;k++) {
//         p[0]=N[k]*N[k];
//         p[1]=N[k+1]*N[k+1];
//         scalar_product+=rho[k]*dz[k]*fe_integraleLGP1xLGP1xLGP1_1D(p, &modes[i][k], &modes[j][k]);
//         }
//       if(verbose==1) printf("W-modes %d %d, scalar product =%lf\n",i,j,scalar_product);
//       }
//     }
//    
//   return 0;
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int normalize_Wmodes(int nlevels, double *z, double **modes, double mask, int nmodes, double *R, double *N, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, k, l;
  int nlayers=nlevels-1;
  double norm;
  double *RR=0, *NN=0;
  double *dz=new double[nlayers];
  
  if(N==0) {
    NN=new double[nlevels];
    for(l=0;l<nlevels;l++) NN[l]=1.0;
    }
  else {
    NN=N;
    }
  if(R==0) {
    RR=new double[nlayers];
    for(k=0;k<nlayers;k++) RR[k]=1.0;
    }
  else {
    RR=R;
    }
    
/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int k=0;k<nlayers;k++) {
    dz[k] = z[k+1]-z[k];
    }
  
//   for(i=0;i<nmodes;i++) {
//     norm=0;
//     double p[2], q[2], r[2], s[2];
//     for (k=1;k<nlevels-1;k++) {
//       norm+=RR[k]*NN[k]*NN[k]*(dz[k-1]+dz[k])*modes[i][k]*modes[i][k]/2.0;
//       }
//     k=0;
//     norm+=RR[k]*NN[k]*NN[k]*dz[k]*modes[i][k]*modes[i][k]/2.0;
//     k=nlevels-1;
//     norm+=RR[k-1]*N[k]*N[k]*dz[k-1]*modes[i][k]*modes[i][k]/2.0;
//     if(norm==0) {
//       if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
//       continue;
//       }
//     norm=sqrt(norm);
//     if(verbose==1) printf("W-mode %d, norm =%lf\n",i,norm);
//     for (k=0;k<nlevels;k++) {
//       modes[i][k]/=norm;
//       }
//     }
    
  for(i=0;i<nmodes;i++) {
    norm=0;
    double p[2], weight=0.0;
    for (k=0;k<nlevels-1;k++) {
//       if(modes[i][k]==mask) continue;
      if(modes[i][k+1]==mask) break;
      p[0]=NN[k]*NN[k];
      p[1]=NN[k+1]*NN[k+1];
      norm+=RR[k]*dz[k]*fe_integraleLGP1xLGP1xLGP1_1D(p, &modes[i][k], &modes[i][k]);
      weight+=dz[k];
      }
    if(norm==0) {
      if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
      continue;
      }
    norm=sqrt(norm/weight);
    if(verbose==1) printf("W-mode %d, norm =%lf\n",i,norm);
    for (k=0;k<nlevels;k++) {
      if(modes[i][k]==mask) continue;
      modes[i][k]/=norm;
      }
    }
  
  if(N==0) delete[] NN;
  if(R==0) delete[] RR;
 
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int normalize_Pmodes(int nlevels, double *z, double **modes, double mask, int nmodes, double *R, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, /*j,*/ k;
  int nlayers=nlevels-1;
  double norm/*,scalar_product*/;
  double *RR=0;
  
  double *dz=new double[nlayers];

  if(R==0) {
    RR=new double[nlayers];
    for(k=0;k<nlayers;k++) RR[k]=1.0;
    }
  else {
    RR=R;
    }
    
/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int k=0;k<nlayers;k++) {
    dz[k] = z[k+1]-z[k];
    }
    
/**----------------------------------------------------------------------------
  normalize eignevectors to unity for integrale norm with weight function*/
  for(i=0;i<nmodes;i++) {
    norm=0;
    double weight=0.0;;
    for (int k=0;k<nlayers;k++) {
      if(modes[i][k]==mask) continue;
      norm+=dz[k]*modes[i][k]*modes[i][k]/RR[k];;
      weight+=dz[k];
      }
    if(norm==0) {
      if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
      continue;
      }
    norm=sqrt(norm/weight);
    if(verbose==1) printf("P-mode %d, norm =%lf\n",i,norm);
    for (int k=0;k<nlayers;k++) {
      if(modes[i][k]==mask) continue;
      modes[i][k]/=norm;
      }
    }
  
//   for(i=0;i<nmodes;i++) {
//     for(j=i;j<nmodes;j++) {
//       scalar_product=0;
//       double p[2], q[2], r[2], s[2];
//       for (int k=0;k<nlevels-1;k++) {
//         scalar_product+=dz[k]*modes[i][k]*modes[j][k]/rho[k];
//         }
//       if(verbose==1) printf("P-modes %d %d, scalar product =%lf\n",i,j,scalar_product);
//       }
//     }
  
  if(R==0) delete[] RR;
 
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int normalize_Umodes(int nlevels, double *z, double **modes, double mask, int nmodes, double *R, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, /*j,*/ k;
  int nlayers=nlevels-1;
  double norm/*,scalar_product*/;
  double *RR=0;
  
  double *dz=new double[nlayers];

  if(R==0) {
    RR=new double[nlayers];
    for(k=0;k<nlayers;k++) RR[k]=1.0;
    }
  else {
    RR=R;
    }
    
/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int k=0;k<nlayers;k++) {
    dz[k] = z[k+1]-z[k];
    }
  
/**----------------------------------------------------------------------------
  normalize eignevectors to unity for integrale norm with weight function*/
  for(i=0;i<nmodes;i++) {
    norm=0;
    double weight=0.0;;
    for (int k=0;k<nlevels-1;k++) {
      if(modes[i][k]==mask) continue;
      norm+=dz[k]*modes[i][k]*modes[i][k]*RR[k];
      weight+=dz[k];
      }
    if(norm==0) {
      if(verbose==1) printf("mode %d, norm =%lf\n",i,norm);
      continue;
      }
    norm=sqrt(norm/weight);
    if(verbose==1) printf("U-mode %d, norm =%lf\n",i,norm);
    for (int k=0;k<nlayers;k++) {
      if(modes[i][k]==mask) continue;
      modes[i][k]/=norm;
      }
    }
  
//   for(i=0;i<nmodes;i++) {
//     for(j=i;j<nmodes;j++) {
//       scalar_product=0;
//       double p[2], q[2], r[2], s[2];
//       for (int k=0;k<nlevels-1;k++) {
//         scalar_product+=dz[k]*modes[i][k]*modes[j][k]*rho[k];
//         }
//       if(verbose==1) printf("U-modes %d %d, scalar product =%lf\n",i,j,scalar_product);
//       }
//     }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_Wmode_v1(int nlevels,double *Z,double *density,double *modes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Solve the eigenvalue problem,
   (equation 34 page 12 of "Developer manuals series Volume IV:
   3D gravity adjustment problem" draft version 24/02/2012
 
*  DSTERF computes all eigenvalues of a symmetric tridiagonal matrix
*  using the Pal-Walker-Kahan variant of the QL or QR algorithm.

*  DSTEVD computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric tridiagonal matrix. If eigenvectors are desired, it
*  uses a divide and conquer algorithm.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.

/// http://www.netlib.org/lapack/double/dstevd.f
/// http://www.netlib.org/lapack/double/dsterf.f
  

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@**/
{
  double Adiag[nlevels],Adown[nlevels-1],N2square[nlevels],dZ_square;
  int LWORK  =1 + 4*nlevels + nlevels*nlevels;
  int LIWORK =3 + 5*nlevels;
  double WORK[LWORK],V[nlevels*nlevels];
  int    IWORK[LIWORK];
  char JOB='V';
  int status=-1;

  for (int i=0;i<nlevels-1;i++) {
/**----------------------------------------------------------------------------
    Z is the level depth*/
    dZ_square=(Z[i+1]-Z[i])*(Z[i+1]-Z[i]);
    Adiag[i]=-2/dZ_square;
    Adown[i]=1/dZ_square;
    }

/**----------------------------------------------------------------------------
  we should apply BC here rigid lid and flat bottom. */
//   dZ_square=(Z[nlevels]-Z[nlevels-1])*(Z[nlevels]-Z[nlevels-1]);
//   Adiag[nlevels]=-2/dZ_square;
  Adiag[0]=0.0;
  Adiag[nlevels-1]=0.0;
  
/**----------------------------------------------------------------------------
  then solve the eigenvalue problem */
  dsterf_(&nlevels,Adiag, Adown, &status);

  printf("After dsterf_ solver: %d  \n",status);
  
  dstevd_(&JOB, &nlevels, Adiag, Adown, V, &nlevels, WORK, IWORK, &LWORK, &LIWORK, &status);

  printf("After dstevd_ solver: %d  \n",status);
  
/**----------------------------------------------------------------------------
  compute N² at element center : N²=(-g/rho dhro /dz) */
  for (int i=0;i<nlevels;i++) {
/**----------------------------------------------------------------------------
    get N² from rho */
    N2square[i]=1.0;
/**----------------------------------------------------------------------------
   then get the velocity modes from the eigen values vector: */
    //modes[i]=sqrt(-1.0*N2square[i]/Adiag[i]);
    modes[i]=1/Adiag[i];
    printf("mode level %d %lf \n",i,modes[i]);
    }
    
return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_Wmode_v2(int nlevels, double *z, double *rho, double g, double *N, double *c, double **modes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Solve the eigenvalue problem,
   (equation 34 page 12 of "Developer manuals series Volume IV:
   3D gravity adjustment problem" draft version 24/02/2012
 
*  DSYGV computes all the eigenvalues, and optionally, the eigenvectors
*  of a real generalized symmetric-definite eigenproblem, of the form
*  A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.
*  Here A and B are assumed to be symmetric and B is also
*  positive definite.

/// http://www.netlib.org/lapack/double/dsygvd.f
  
*  Arguments
*  =========
*
*  ITYPE   (input) INTEGER
*          Specifies the problem type to be solved:
*          = 1:  A*x = (lambda)*B*x
*          = 2:  A*B*x = (lambda)*x
*          = 3:  B*A*x = (lambda)*x
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangles of A and B are stored;
*          = 'L':  Lower triangles of A and B are stored.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          matrix Z of eigenvectors.  The eigenvectors are normalized
*          as follows:
*          if ITYPE = 1 or 2, Z**T*B*Z = I;
*          if ITYPE = 3, Z**T*inv(B)*Z = I.
*          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
*          or the lower triangle (if UPLO='L') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
*          On entry, the symmetric matrix B.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of B contains the
*          upper triangular part of the matrix B.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of B contains
*          the lower triangular part of the matrix B.
*
*          On exit, if INFO <= N, the part of B containing the matrix is
*          overwritten by the triangular factor U or L from the Cholesky
*          factorization B = U**T*U or B = L*L**T.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (max(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If N <= 1,               LWORK >= 1.
*          If JOBZ = 'N' and N > 1, LWORK >= 2*N+1.
*          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK and IWORK
*          arrays, returns these values as the first entries of the WORK
*          and IWORK arrays, and no error message related to LWORK or
*          LIWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (max(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If N <= 1,                LIWORK >= 1.
*          If JOBZ  = 'N' and N > 1, LIWORK >= 1.
*          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK and IWORK arrays, and no error message related to
*          LWORK or LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  DPOTRF or DSYEVD returned an error code:
*             <= N:  if INFO = i and JOBZ = 'N', then the algorithm
*                    failed to converge; i off-diagonal elements of an
*                    intermediate tridiagonal form did not converge to
*                    zero;
*                    if INFO = i and JOBZ = 'V', then the algorithm
*                    failed to compute an eigenvalue while working on
*                    the submatrix lying in rows and columns INFO/(N+1)
*                    through mod(INFO,N+1);
*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
*                    minor of order i of B is not positive definite.
*                    The factorization of B could not be completed and
*                    no eigenvalues or eigenvectors were computed.
*
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ **/
{
  int status=-1;
  int k,l,row,col;
  int nlayers=nlevels-1;
  double *dZ,*dz,*Z,*N2;
  int dim=nlevels-1;
  
  int LWORK  =1 + 6*dim + 2*dim*dim;
  int LIWORK =3 + 5*dim;
  int ITYPE=1;
  char JOB='V', UPLO='U';

//#define STATICS

#ifdef STATICS
   double WORK[LWORK],W[dim];
   int IWORK[LIWORK];
   double A[dim*dim],B[dim*dim];
#else
  double *WORK,*W;
  int    *IWORK;
  double *A,*B;
  
  WORK =     new double [LWORK];
  IWORK =    new int    [LIWORK];
  W =        new double [dim];
  A =        new double [dim*dim];
  B =        new double [dim*dim];
#endif
  dz=new double[nlayers];
  Z =new double[nlayers];
  dZ=new double[nlevels];
  N2=new double[nlevels];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int i=0;i<nlayers;i++) {
    dz[i] = z[i+1]-z[i];
    Z[i]  = 0.5*(z[i+1]+z[i]);
    }

/**----------------------------------------------------------------------------
  Z is the layer depth; dZ is "level" thickness*/
  for (int i=1;i<nlevels-1;i++) {
    dZ[i]=Z[i]-Z[i-1];
    }
  dZ[0]=1.e+35;
  dZ[nlevels-1]=1.e+35;
  
/**----------------------------------------------------------------------------
  compute N² at levels : N²=(-g/rho dhro /dz) */
  for (l=1;l<nlevels-1;l++) {
    k=l;
    double r=0.5*(rho[k]+rho[k-1]);
    N2[l]=-g*(rho[k]-rho[k-1])/dZ[l]/r;
    if(N2[l]<1.0e-8) {
//       printf("unstable stratification %lf, abort...\n",N2[l]);
//       return(-1);
      N2[l]=1.0e-8;
      }
    }

  N2[0]=N2[1];
  N2[nlevels-1]=N2[nlevels-2];
  
//   A=new double[dim*dim];
//   B=new double[dim*dim];
  
  for(k=0;k<dim*dim;k++) A[k]=0.0;
  for(k=0;k<dim*dim;k++) B[k]=0.0;

/**----------------------------------------------------------------------------
  skip the 2 first levels (bottom level and next) and last level (surface level)*/
  for (int l=2;l<nlevels-1;l++) {
    k=l;
    row=l-1;
    col=row-1;
/**----------------------------------------------------------------------------
    A matrix */
    A[dim*col+row] = -1/dz[k-1];
    col=row;
    A[dim*col+row] =  1/dz[k-1]+1/dz[k];
    col=row+1;
    A[dim*col+row] = -1/dz[k];
/**----------------------------------------------------------------------------
    B matrix */
    col=row;
    B[dim*col+row] = N2[l]*dZ[l];
    }

/**----------------------------------------------------------------------------
  first level above bottom */
  l=1;
  k=l;
  row=l-1;
  col=row;
  A[dim*col+row] =  1/dz[k-1]+1/dz[k];
  col=row+1;
  A[dim*col+row] = -1/dz[k];
  col=row;
  B[dim*col+row] = N2[l]*dZ[l];
  
/**----------------------------------------------------------------------------
  top level */
  l=nlevels-1;
  k=l;
  row=l-1;
  col=row-1;
  col=row-1;
  A[dim*col+row] = -1/dz[k-1];
  col=row;
  A[dim*col+row] =  1/dz[k-1];
  col=row;
  B[dim*col+row] = g;
  
/**----------------------------------------------------------------------------
  then solve the eigenvalue problem */
  dsygvd_(&ITYPE, &JOB, &UPLO, &dim, A, &dim, B, &dim, W, WORK, &LWORK, IWORK, &LIWORK, &status);

  if(status!=0) {
//    printf("After dsygvd solver: status=%d  dim=%d\n",status, dim);
    goto leave;
    }
  
  for (int i=0;i<nlevels-1;i++) {
/**----------------------------------------------------------------------------
    then get the velocity modes from the eigen values vector: */
    if(W[i]==0) {
      printf("zero eigenvalue %dth : %lf \n",i,W[i]);
      }
    c[i]=1./sqrt(W[i]);
//    STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
    for (int j=0;j<nlevels-1;j++) {
      modes[i][j+1]=A[dim*i+j];
      }
    modes[i][0]=0;
    }
  for (int i=0;i<nlevels;i++) {
    N[i]=sqrt(N2[i]);
    }

leave:
  delete[] dz;
  delete[] dZ;
  delete[] Z;
  delete[] N2;
  
#ifdef STATICS
#else
  delete[] WORK;
  delete[] IWORK;
  delete[] W;
  delete[] A;
  delete[] B;
#endif
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_Wmode_v3(int nlevels,double *z, double *rho, double g, double *N, double *c, double **modes, bool rigid_lid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
*
      SUBROUTINE DGEEV( JOBVL, JOBVR, N, A, LDA, WR, WI, VL, LDVL, VR,
     $                  LDVR, WORK, LWORK, INFO )
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WI( * ), WORK( * ), WR( * )
*     ..
*
*  Purpose
*  =======
*
*  DGEEV computes for an N-by-N real nonsymmetric matrix A, the
*  eigenvalues and, optionally, the left and/or right eigenvectors.
*
*  The right eigenvector v(j) of A satisfies
*                   A * v(j) = lambda(j) * v(j)
*  where lambda(j) is its eigenvalue.
*  The left eigenvector u(j) of A satisfies
*                u(j)**T * A = lambda(j) * u(j)**T
*  where u(j)**T denotes the transpose of u(j).
*
*  The computed eigenvectors are normalized to have Euclidean norm
*  equal to 1 and largest component real.
*
*  Arguments
*  =========
*
*  JOBVL   (input) CHARACTER*1
*          = 'N': left eigenvectors of A are not computed;
*          = 'V': left eigenvectors of A are computed.
*
*  JOBVR   (input) CHARACTER*1
*          = 'N': right eigenvectors of A are not computed;
*          = 'V': right eigenvectors of A are computed.
*
*  N       (input) INTEGER
*          The order of the matrix A. N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N matrix A.
*          On exit, A has been overwritten.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  WR      (output) DOUBLE PRECISION array, dimension (N)
*  WI      (output) DOUBLE PRECISION array, dimension (N)
*          WR and WI contain the real and imaginary parts,
*          respectively, of the computed eigenvalues.  Complex
*          conjugate pairs of eigenvalues appear consecutively
*          with the eigenvalue having the positive imaginary part
*          first.
*
*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
*          If JOBVL = 'V', the left eigenvectors u(j) are stored one
*          after another in the columns of VL, in the same order
*          as their eigenvalues.
*          If JOBVL = 'N', VL is not referenced.
*          If the j-th eigenvalue is real, then u(j) = VL(:,j),
*          the j-th column of VL.
*          If the j-th and (j+1)-st eigenvalues form a complex
*          conjugate pair, then u(j) = VL(:,j) + i*VL(:,j+1) and
*          u(j+1) = VL(:,j) - i*VL(:,j+1).
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.  LDVL >= 1; if
*          JOBVL = 'V', LDVL >= N.
*
*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
*          If JOBVR = 'V', the right eigenvectors v(j) are stored one
*          after another in the columns of VR, in the same order
*          as their eigenvalues.
*          If JOBVR = 'N', VR is not referenced.
*          If the j-th eigenvalue is real, then v(j) = VR(:,j),
*          the j-th column of VR.
*          If the j-th and (j+1)-st eigenvalues form a complex
*          conjugate pair, then v(j) = VR(:,j) + i*VR(:,j+1) and
*          v(j+1) = VR(:,j) - i*VR(:,j+1).
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.  LDVR >= 1; if
*          JOBVR = 'V', LDVR >= N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (max(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,3*N), and
*          if JOBVL = 'V' or JOBVR = 'V', LWORK >= 4*N.  For good
*          performance, LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, the QR algorithm failed to compute all the
*                eigenvalues, and no eigenvectors have been computed;
*                elements i+1:N of WR and WI contain eigenvalues which
*                have converged.
*
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@**/
{
  int status=-1;
  int k,l,row,col;
  int nlayers=nlevels-1;
  double *dZ,*dz,*Z,*N2;
  int dim=nlevels-1,neq=dim;
  
  double unsorted[dim], A[dim*dim], B[dim*dim];
  int LWORK  =1 + 4*dim;
  double WORK[LWORK],WR[dim],WI[dim],VR[dim*dim],VL[dim*dim];
  
  dz=new double[nlayers];
  Z =new double[nlayers];
  dZ=new double[nlevels];
  N2=new double[nlevels];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int i=0;i<nlayers;i++) {
    dz[i] = z[i+1]-z[i];
    Z[i]  = 0.5*(z[i+1]+z[i]);
    }

/**----------------------------------------------------------------------------
  Z is the layer depth; dZ is "level" thickness*/
  for (int i=1;i<nlevels-1;i++) {
    dZ[i]=Z[i]-Z[i-1];
    }
  dZ[0]=1.e+35;
  dZ[nlevels-1]=1.e+35;
  
/**----------------------------------------------------------------------------
  compute N² at levels : N²=(-g/rho dhro /dz) */
  for (l=1;l<nlevels-1;l++) {
    k=l;
    double r=0.5*(rho[k]+rho[k-1]);
    N2[l]=-g*(rho[k]-rho[k-1])/dZ[l]/r;
//    N2[l]=2.e-03*2.e-03;
    N2[l]=max(N2[l], 1.e-10);
    if(N2[l]<0) {
      printf("compute_Wmode_v3: unstable stratification, abort...\n");
      return(-1);
      }
    }

  N2[0]=N2[1];
  N2[nlevels-1]=N2[nlevels-2];
  
  for(k=0;k<dim*dim;k++) A[k]=0.0;
  for(k=0;k<dim*dim;k++) B[k]=0.0;

/**----------------------------------------------------------------------------
  skip the 2 first levels (bottom level and next) and last level (surface level)*/
  for (int l=2;l<nlevels-1;l++) {
    k=l;
    row=l-1;
/**----------------------------------------------------------------------------
    A matrix */
    col=row-1;
    A[dim*col+row] = -1.0/dz[k-1];
    col=row;
    A[dim*col+row] =  1.0/dz[k-1]+1/dz[k];
    col=row+1;
    A[dim*col+row] = -1.0/dz[k];
/**----------------------------------------------------------------------------
    B matrix */
    col=row;
    B[dim*col+row] = N2[l]*dZ[l];
    }

/**----------------------------------------------------------------------------
  first level above bottom */
  l=1;
  k=l;
  row=l-1;
  col=row;
  A[dim*col+row] =  1/dz[k-1]+1/dz[k];
  col=row+1;
  A[dim*col+row] = -1/dz[k];
  col=row;
  B[dim*col+row] = N2[l]*dZ[l];
  
/**----------------------------------------------------------------------------
  top level */
  l=nlevels-1;
  k=l;
  row=l-1;
  
  rigid_lid=false;
  
  if(rigid_lid) {
    neq--;
    }
  else {
    col=row-1;
    A[dim*col+row] = -1/dz[k-1];
    col=row;
    A[dim*col+row] =  1/dz[k-1];
    col=row;
    B[dim*col+row] = g;
    }
  
/**----------------------------------------------------------------------------
  perform inv(B)*A */
  for (int row=0;row<nlevels-1;row++) {
/**----------------------------------------------------------------------------
    B matrix */
    col=row;
    double r=B[dim*col+row];
    for (int col=0;col<nlevels-1;col++) {
      A[dim*col+row]/=r;
      }
    }
    
/**----------------------------------------------------------------------------
  then solve the eigenvalue problem; imaginary part issue to be discussed */
  char JOBVL='N',JOBVR='V';
  dgeev_( &JOBVL, &JOBVR, &neq, A, &dim, WR, WI, VL, &dim, VR, &dim, WORK, &LWORK, &status );
  
  if(verbose==1) STDOUT_BASE_LINE("After dgeev solver: %d  \n",status);
  
  if(status!=0) return(-1);
  
  for (int k=0;k<nlevels-1;k++) {
    for (int j=0;j<nlevels-1;j++) {
      modes[k][j]=0;
      }
    }
  
/**----------------------------------------------------------------------------
  Enforce decreasing order (not granted by solver) */
  for (int k=0;k<neq;k++) {
/**----------------------------------------------------------------------------
    get the velocity modes from the eigen values vector */
    unsorted[k]=1./sqrt(WR[k]);
    }
  size_t *order=sort(unsorted,neq);
  
  for (int i=0;i<neq;i++) {
    int k=order[neq-1-i];
    c[i]=1./sqrt(WR[k]);
    if(isfinite(c[i])==0) {
      STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
      }
    if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
    modes[i][dim-1]=0;
/**----------------------------------------------------------------------------
    eigenvector from second level to top level */
    for (int j=0;j<neq;j++) {
      if(isfinite(VR[dim*k+j])==0) {
        STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
        }
      modes[i][j+1]=VR[dim*k+j];
      }
/**----------------------------------------------------------------------------
    and set bottom value to zero (flat bottom boundary condition)*/
    modes[i][0]=0;
    }
    
  for (int i=0;i<nlevels;i++) {
    N[i]=sqrt(N2[i]);
    }
    
  delete[] dz;
  delete[] dZ;
  delete[] Z;
  delete[] N2;
  
  delete[] order;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_Wmode_v4(int nlevels,double *z, double *rho, double g, double *N, double *c, double **modes, bool rigid_lid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  int k,l,row,col;
  int nlayers=nlevels-1;
  double *dZ,*dz,*Z,*N2;
  int dim=nlevels-1,neq=dim;
  
  double unsorted[dim], A[dim*dim], B[dim*dim];
  int LWORK  =1 + 4*dim;
  double WORK[LWORK],WR[dim],WI[dim],VR[dim*dim],VL[dim*dim];
  
  dz=new double[nlayers];
  Z =new double[nlayers];
  dZ=new double[nlevels];
  N2=new double[nlevels];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int i=0;i<nlayers;i++) {
    dz[i] = z[i+1]-z[i];
    Z[i]  = 0.5*(z[i+1]+z[i]);
    }

/**----------------------------------------------------------------------------
  Z is the layer depth; dZ is "level" thickness*/
  for (int i=1;i<nlevels-1;i++) {
    dZ[i]=Z[i]-Z[i-1];
    }
  dZ[0]=1.e+35;
  dZ[nlevels-1]=1.e+35;
  
/**----------------------------------------------------------------------------
  compute N² at levels : N²=(-g/rho dhro /dz) */
  for (l=1;l<nlevels-1;l++) {
    k=l;
    double r=0.5*(rho[k]+rho[k-1]);
    N2[l]=-g*(rho[k]-rho[k-1])/dZ[l]/r;
    N2[l]=max(N2[l], 1.e-15);
    if(N2[l]<0) {
      printf("compute_Wmode_v4: unstable stratification, abort...\n");
      return(-1);
      }
    }

  N2[0]=N2[1];
  N2[nlevels-1]=N2[nlevels-2];
  
  for(k=0;k<dim*dim;k++) A[k]=0.0;
  for(k=0;k<dim*dim;k++) B[k]=0.0;

/**----------------------------------------------------------------------------
  skip the 2 first levels (bottom level and next) and last level (surface level)*/
  for (int l=2;l<nlevels-1;l++) {
    k=l;
    row=l-1;
/**----------------------------------------------------------------------------
    A matrix */
    col=row-1;
    A[dim*col+row] = -1.0*rho[k-1]/dz[k-1];
    col=row;
    A[dim*col+row] =  1.0*rho[k-1]/dz[k-1]+1.0*rho[k]/dz[k];
    col=row+1;
    A[dim*col+row] = -1.0*rho[k]/dz[k];
/**----------------------------------------------------------------------------
    B matrix */
    double r=0.5*(rho[k]+rho[k-1]);
    col=row;
    B[dim*col+row] = r*N2[l]*dZ[l];
    }

/**----------------------------------------------------------------------------
  first level above bottom */
  l=1;
  k=l;
  row=l-1;
  col=row;
  A[dim*col+row] =  rho[k-1]/dz[k-1]+rho[k]/dz[k];
  col=row+1;
  A[dim*col+row] = -rho[k]/dz[k];
  col=row;
  double r=0.5*(rho[k]+rho[k-1]);
  B[dim*col+row] = r*N2[l]*dZ[l];
  
/**----------------------------------------------------------------------------
  top level */

  rigid_lid=false;
  
  l=nlevels-1;
  k=l;
  row=l-1;
  
  if(rigid_lid) {
    neq--;
    }
  else {
    col=row-1;
    A[dim*col+row] = -1/dz[k-1];
    col=row;
    A[dim*col+row] =  1/dz[k-1];
    col=row;
    B[dim*col+row] = g;
    }
  
/**----------------------------------------------------------------------------
  perform inv(B)*A */
  for (int row=0;row<nlevels-1;row++) {
/**----------------------------------------------------------------------------
    B matrix */
    col=row;
    double r=B[dim*col+row];
    for (int col=0;col<nlevels-1;col++) {
      A[dim*col+row]/=r;
      }
    }
    
/**----------------------------------------------------------------------------
  then solve the eigenvalue problem; imaginary part issue to be discussed */
  char JOBVL='N',JOBVR='V';
  dgeev_( &JOBVL, &JOBVR, &neq, A, &dim, WR, WI, VL, &dim, VR, &dim, WORK, &LWORK, &status );
  
  if(verbose==1) STDOUT_BASE_LINE("After dgeev solver: %d  \n",status);
  
  if(status!=0) return(-1);
  
  for (int k=0;k<nlevels-1;k++) {
    for (int j=0;j<nlevels-1;j++) {
      modes[k][j]=0;
      }
    }
  
/**----------------------------------------------------------------------------
  Enforce decreasing order (not granted by solver) */
  for (int k=0;k<neq;k++) {
/**----------------------------------------------------------------------------
    get the velocity modes from the eigen values vector */
    unsorted[k]=1./sqrt(WR[k]);
    }
  size_t *order=sort(unsorted,neq);
  
  for (int i=0;i<neq;i++) {
    int k=order[neq-1-i];
    c[i]=1./sqrt(WR[k]);
    if(isfinite(c[i])==0) {
      STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
      }
    if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
    modes[i][dim-1]=0;
/**----------------------------------------------------------------------------
    eigenvector from second level to top level */
    for (int j=0;j<neq;j++) {
      if(isfinite(VR[dim*k+j])==0) {
        STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
        }
      modes[i][j+1]=VR[dim*k+j];
      }
/**----------------------------------------------------------------------------
    and set bottom value to zero (flat bottom boundary condition)*/
    modes[i][0]=0;
    }
    
  for (int i=0;i<nlevels;i++) {
    N[i]=sqrt(N2[i]);
    }
    
  delete[] dz;
  delete[] dZ;
  delete[] Z;
  delete[] N2;
  
  delete[] order;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_Wmode_v5(int nlevels, double *z, double *rho, double g, double *N, double *c, double **modes, int & nmodes, bool rigid_lid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*SUBROUTINE DGGEV( JOBVL, JOBVR, N, A, LDA, B, LDB, ALPHAR, ALPHAI,
     $                  BETA, VL, LDVL, VR, LDVR, WORK, LWORK, INFO )
*
*  -- LAPACK driver routine (version 3.2) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2006
*
*     .. Scalar Arguments ..
      CHARACTER          JOBVL, JOBVR
      INTEGER            INFO, LDA, LDB, LDVL, LDVR, LWORK, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), ALPHAI( * ), ALPHAR( * ),
     $                   B( LDB, * ), BETA( * ), VL( LDVL, * ),
     $                   VR( LDVR, * ), WORK( * )
*     ..
*
*  Purpose
*  =======
*
*  DGGEV computes for a pair of N-by-N real nonsymmetric matrices (A,B)
*  the generalized eigenvalues, and optionally, the left and/or right
*  generalized eigenvectors.
*
*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar
*  lambda or a ratio alpha/beta = lambda, such that A - lambda*B is
*  singular. It is usually represented as the pair (alpha,beta), as
*  there is a reasonable interpretation for beta=0, and even for both
*  being zero.
*
*  The right eigenvector v(j) corresponding to the eigenvalue lambda(j)
*  of (A,B) satisfies
*
*                   A * v(j) = lambda(j) * B * v(j).
*
*  The left eigenvector u(j) corresponding to the eigenvalue lambda(j)
*  of (A,B) satisfies
*
*                   u(j)**H * A  = lambda(j) * u(j)**H * B .
*
*  where u(j)**H is the conjugate-transpose of u(j).
*
*
*  Arguments
*  =========
*
*  JOBVL   (input) CHARACTER*1
*          = 'N':  do not compute the left generalized eigenvectors;
*          = 'V':  compute the left generalized eigenvectors.
*
*  JOBVR   (input) CHARACTER*1
*          = 'N':  do not compute the right generalized eigenvectors;
*          = 'V':  compute the right generalized eigenvectors.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, VL, and VR.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the matrix A in the pair (A,B).
*          On exit, A has been overwritten.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
*          On entry, the matrix B in the pair (A,B).
*          On exit, B has been overwritten.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  LDB >= max(1,N).
*
*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
*          be the generalized eigenvalues.  If ALPHAI(j) is zero, then
*          the j-th eigenvalue is real; if positive, then the j-th and
*          (j+1)-st eigenvalues are a complex conjugate pair, with
*          ALPHAI(j+1) negative.
*
*          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
*          may easily over- or underflow, and BETA(j) may even be zero.
*          Thus, the user should avoid naively computing the ratio
*          alpha/beta.  However, ALPHAR and ALPHAI will be always less
*          than and usually comparable with norm(A) in magnitude, and
*          BETA always less than and usually comparable with norm(B).
*
*  VL      (output) DOUBLE PRECISION array, dimension (LDVL,N)
*          If JOBVL = 'V', the left eigenvectors u(j) are stored one
*          after another in the columns of VL, in the same order as
*          their eigenvalues. If the j-th eigenvalue is real, then
*          u(j) = VL(:,j), the j-th column of VL. If the j-th and
*          (j+1)-th eigenvalues form a complex conjugate pair, then
*          u(j) = VL(:,j)+i*VL(:,j+1) and u(j+1) = VL(:,j)-i*VL(:,j+1).
*          Each eigenvector is scaled so the largest component has
*          abs(real part)+abs(imag. part)=1.
*          Not referenced if JOBVL = 'N'.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the matrix VL. LDVL >= 1, and
*          if JOBVL = 'V', LDVL >= N.
*
*  VR      (output) DOUBLE PRECISION array, dimension (LDVR,N)
*          If JOBVR = 'V', the right eigenvectors v(j) are stored one
*          after another in the columns of VR, in the same order as
*          their eigenvalues. If the j-th eigenvalue is real, then
*          v(j) = VR(:,j), the j-th column of VR. If the j-th and
*          (j+1)-th eigenvalues form a complex conjugate pair, then
*          v(j) = VR(:,j)+i*VR(:,j+1) and v(j+1) = VR(:,j)-i*VR(:,j+1).
*          Each eigenvector is scaled so the largest component has
*          abs(real part)+abs(imag. part)=1.
*          Not referenced if JOBVR = 'N'.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the matrix VR. LDVR >= 1, and
*          if JOBVR = 'V', LDVR >= N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (max(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,8*N).
*          For good performance, LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          = 1,...,N:
*                The QZ iteration failed.  No eigenvectors have been
*                calculated, but ALPHAR(j), ALPHAI(j), and BETA(j)
*                should be correct for j=INFO+1,...,N.
*          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
*                =N+2: error return from DTGEVC.
*
*  =====================================================================*/
{
  int status=-1;
  int k,l,row,col;
  int nlayers=nlevels-1;
  double *dZ,*dz,*Z,*N2;
  int dim=nlevels-1,neq=dim;
  
  double A[dim*dim],B[dim*dim];
  int LWORK  =1 + 8*dim;
  double WORK[LWORK],ALPHAR[dim],ALPHAI[dim],BETA[dim],VR[dim*dim],VL[dim*dim];
  double unsorted[dim];
 
  dz=new double[nlayers];
  Z =new double[nlayers];
  dZ=new double[nlevels];
  N2=new double[nlevels];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int i=0;i<nlayers;i++) {
    dz[i] = z[i+1]-z[i];
    Z[i]  = 0.5*(z[i+1]+z[i]);
    }

/**----------------------------------------------------------------------------
  Z is the layer depth; dZ is "level" thickness*/
  for (int i=1;i<nlevels-1;i++) {
    dZ[i]=Z[i]-Z[i-1];
    }
  dZ[0]=1.e+35;
  dZ[nlevels-1]=1.e+35;
  
/**----------------------------------------------------------------------------
  compute N² at levels : N²=(-g/rho dhro /dz) */
  for (l=1;l<nlevels-1;l++) {
    k=l;
    double r=0.5*(rho[k]+rho[k-1]);
    N2[l]=-g * (rho[k]-rho[k-1])/dZ[l]/r;
    N2[l]=max(N2[l], 1.e-40);
    if(N2[l]<0) {
      printf("compute_Wmode_v5: unstable stratification, abort...\n");
      return(-1);
      }
    }

  N2[0]=N2[1];
  N2[nlevels-1]=N2[nlevels-2];
  
  for(k=0;k<dim*dim;k++) A[k]=0.0;
  for(k=0;k<dim*dim;k++) B[k]=0.0;

/**----------------------------------------------------------------------------
  skip the 2 first levels (bottom level and next) and last level (surface level)*/
  for (int l=2;l<nlevels-1;l++) {
    k=l;
    row=l-1;
/**----------------------------------------------------------------------------
    A matrix */
    col=row-1;
    A[dim*col+row] = -rho[k-1]/dz[k-1];
    col=row;
    A[dim*col+row] =  rho[k-1]/dz[k-1]+rho[k]/dz[k];
    col=row+1;
    A[dim*col+row] = -rho[k]/dz[k];
/**----------------------------------------------------------------------------
    B matrix */
    double r=0.5*(rho[k]+rho[k-1]);
    col=row;
    B[dim*col+row] = r*N2[l]*dZ[l];
    }

/**----------------------------------------------------------------------------
  first level above bottom */
  l=1;
  k=l;
  row=l-1;
  col=row;
  A[dim*col+row] =  rho[k-1]/dz[k-1]+rho[k]/dz[k];
  col=row+1;
  A[dim*col+row] = -rho[k]/dz[k];
  col=row;
  double r=0.5*(rho[k]+rho[k-1]);
  B[dim*col+row] = r*N2[l]*dZ[l];
  
/**----------------------------------------------------------------------------
  top level */
  l=nlevels-1;
  k=l;
  row=l-1;
  
  if(rigid_lid) {
    neq--;
    }
  else {
    col=row-1;
    A[dim*col+row] = -1/dz[k-1];
    col=row;
    A[dim*col+row] =  1/dz[k-1];
    col=row;
    B[dim*col+row] = g;
    }
      
/**----------------------------------------------------------------------------
  then solve the eigenvalue problem; imaginary part issue to be discussed */
  char JOBVL='N',JOBVR='V';
  dggev_( &JOBVL, &JOBVR, &neq, A, &dim, B, &dim, ALPHAR, ALPHAI, BETA, VL, &dim, VR, &dim, WORK, &LWORK, &status );
  
  if(verbose==1) STDOUT_BASE_LINE("After dgeev solver: %d  \n",status);
  
  if(status!=0) return(-1);
  
  for (int k=0;k<nlevels-1;k++) {
    for (int j=0;j<nlevels-1;j++) {
      modes[k][j]=0;
      }
    }
  
  nmodes=neq;
  for (int k=0;k<neq;k++) {
    if(fabs(BETA[k]) > 1.e-10) ALPHAR[k]/=BETA[k];
    else {
      ALPHAR[k]=1.e+20;
      nmodes--;
      }
    }
  
/**----------------------------------------------------------------------------
  Enforce decreasing order (not granted by solver) */
  for (int k=0;k<neq;k++) {
/**----------------------------------------------------------------------------
    get the velocity modes from the eigen values vector */
    unsorted[k]=1./sqrt(ALPHAR[k]);
    }
  size_t *order=sort(unsorted,neq);
  
  for (int i=0;i<neq;i++) {
    int k=order[neq-1-i];
    c[i]=1./sqrt(ALPHAR[k]);
    if(ALPHAR[k]==1.e+20) {
      for (int j=0;j<neq+1;j++) {
        modes[i][j]=0;
        }
      continue;
      }
    if(isfinite(c[i])==0) {
      STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
      }
    if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
    modes[i][dim-1]=0;
/**----------------------------------------------------------------------------
    eigenvector from second level to top level */
    for (int j=0;j<neq;j++) {
      if(isfinite(VR[dim*k+j])==0) {
        if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
        }
      modes[i][j+1]=VR[dim*k+j];
      }
/**----------------------------------------------------------------------------
    and set bottom value to zero (flat bottom boundary condition)*/
    modes[i][0]=0;
    }
    
  for (int i=0;i<nlevels;i++) {
    N[i]=sqrt(N2[i]);
    }
    
  delete[] dz;
  delete[] dZ;
  delete[] Z;
  delete[] N2;
  
  delete[] order;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_Wmode_v6(int nlevels, double *z, double *T, double *S, double *rho, double g, double *N, double *c, double **modes, int & nmodes, bool rigid_lid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  int k,l,row,col;
  const int nlayers=nlevels-1;
  double *dZ,*dz,*N2;
  int dim=nlevels,neq=dim;
  
  double A[dim*dim],B[dim*dim];
  int LWORK  =1 + 8*dim;
  double WORK[LWORK],ALPHAR[dim],ALPHAI[dim],BETA[dim],VR[dim*dim],VL[dim*dim];
  double unsorted[dim];
  
  double z_layer[nlayers];
 
  dz=new double[nlayers];
  dZ=new double[nlevels];
  N2=new double[nlevels];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int i=0;i<nlayers;i++) {
    dz[i] = z[i+1]-z[i];
    z_layer[i]  = 0.5*(z[i+1]+z[i]);
    }

/**----------------------------------------------------------------------------
  z_layer is the layer depth; dZ is "level" thickness*/
  for (int i=1;i<nlevels-1;i++) {
    dZ[i]=z_layer[i]-z_layer[i-1];
    }
  dZ[0]=1.e+35;
  dZ[nlevels-1]=1.e+35;
  
/**----------------------------------------------------------------------------
  compute N² at levels : N²=(-g/rho dhro /dz) */
  for (l=1;l<nlevels-1;l++) {
    k=l;
    if(T==0) {
      double r=0.5*(rho[k]+rho[k-1]);
      N2[l]=-g * (rho[k]-rho[k-1])/dZ[l]/r;
      N2[l]=max(N2[l], 1.e-40);
      }
    else {
      double r1=water_density(T[k-1], S[k-1], -z[l], 1);
      double r2=water_density(T[k],   S[k],   -z[l], 1);
      double rlocal=0.5*(r1+r2);
      N2[l]=-g*(r2-r1)/dZ[l]/rlocal;
      N2[l]=max(N2[l], 0.0);
      }
//     if(N2[l]<0) {
//       printf("compute_Wmode_v6: unstable stratification, abort...\n");
//       return(-1);
//       }
    }

  N2[0]=N2[1];
  N2[nlevels-1]=N2[nlevels-2];
  
  for(k=0;k<dim*dim;k++) A[k]=0.0;
  for(k=0;k<dim*dim;k++) B[k]=0.0;

/**----------------------------------------------------------------------------
  skip the first level (bottom level) and last level (surface level)*/
  for (int l=1;l<nlevels-1;l++) {
    k=l;
    row=l;
    double r1,r2;
    if(T==0) {
      r1=rho[k-1];
      r2=rho[k];
      }
    else {
      r1=water_density(T[k-1], S[k-1], -z[l], 1);
      r2=water_density(T[k],   S[k],   -z[l], 1);
      }
/**----------------------------------------------------------------------------
    A matrix */
    col=row-1;
    A[dim*col+row] = -r1/dz[k-1];
    col=row;
    A[dim*col+row] =  r1/dz[k-1]+r2/dz[k];
    col=row+1;
    A[dim*col+row] = -r2/dz[k];
/**----------------------------------------------------------------------------
    B matrix */
    double r=0.5*(r2+r1);
    col=row;
    B[dim*col+row] = r*N2[l]*dZ[l];
    }

/**----------------------------------------------------------------------------
  bottom level*/
  l=0;
  k=0;
  row=l;
  col=row;
  A[dim*col+row] =  1.0;
  col=row;
  B[dim*col+row] =  0.0;
  
/**----------------------------------------------------------------------------
  top level */
  l=nlevels-1;
  k=l;
  row=l;
  
  if(rigid_lid) {
    neq--;
    }
  else {
    col=row-1;
    A[dim*col+row] = -1/dz[k-1];
    col=row;
    A[dim*col+row] =  1/dz[k-1];
    col=row;
    B[dim*col+row] = g;
    }
      
/**----------------------------------------------------------------------------
  then solve the eigenvalue problem; imaginary part issue to be discussed */
  char JOBVL='N',JOBVR='V';
  dggev_( &JOBVL, &JOBVR, &neq, A, &dim, B, &dim, ALPHAR, ALPHAI, BETA, VL, &dim, VR, &dim, WORK, &LWORK, &status );
  
  if(verbose==1) STDOUT_BASE_LINE("After dgeev solver: %d  \n",status);
  
  if(status!=0) return(-1);
  
  for (int k=0;k<nlevels;k++) {
    for (int j=0;j<nlevels;j++) {
      modes[k][j]=0;
      }
    }
  
  nmodes=neq;
  for (int k=0;k<neq;k++) {
    if(fabs(BETA[k]) > 1.e-10) ALPHAR[k]/=BETA[k];
    else {
      ALPHAR[k]=1.e+20;
      nmodes--;
      }
    }
  
/**----------------------------------------------------------------------------
  Enforce decreasing order (not granted by solver) */
  for (int k=0;k<neq;k++) {
/**----------------------------------------------------------------------------
    get the velocity modes from the eigen values vector */
    unsorted[k]=1./sqrt(ALPHAR[k]);
    }
  size_t *order=sort(unsorted,neq);
  
  for (int i=0;i<neq;i++) {
    int k=order[neq-1-i];
    c[i]=1./sqrt(ALPHAR[k]);
    if(ALPHAR[k]==1.e+20) {
      for (int j=0;j<neq;j++) {
        modes[i][j]=0;
        }
      continue;
      }
    if( not isnormal(c[i]) ) {
      STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
      }
    if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
    modes[i][dim-1]=0;
/**----------------------------------------------------------------------------
    eigenvector from second level to top level */
    for (int j=0;j<neq;j++) {
      if(isfinite(VR[dim*k+j])==0) {
        if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
        }
      modes[i][j]=VR[dim*k+j];
      }
    }
    
  for (int i=0;i<nlevels;i++) {
    N[i]=sqrt(N2[i]);
    }
    
  delete[] dz;
  delete[] dZ;
  delete[] N2;
  
  delete[] order;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_Pmode(int nlevels,double *z, double *rho, double g, double *N, double *c, double **modes, int & nmodes, bool rigid_lid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  int k,l,row,col;
  int nlayers=nlevels-1;
  double *dZ,*dz,*Z,*N2;
  int dim=nlevels-1,neq=dim;
  
  double A[dim*dim],B[dim*dim];
  int LWORK  =1 + 8*dim;
  double WORK[LWORK],ALPHAR[dim],ALPHAI[dim],BETA[dim],VR[dim*dim],VL[dim*dim];
  double unsorted[dim];
  
  bool debug=false;
  
  dz=new double[nlayers];
  Z =new double[nlayers];
  dZ=new double[nlevels];
  N2=new double[nlevels];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int i=0;i<nlayers;i++) {
    dz[i] = z[i+1]-z[i];
    Z[i]  = 0.5*(z[i+1]+z[i]);
    }

/**----------------------------------------------------------------------------
  Z is the layer depth; dZ is "level" thickness*/
  for (int i=1;i<nlevels-1;i++) {
    dZ[i]=Z[i]-Z[i-1];
    }
  dZ[0]=1.e+35;
  dZ[nlevels-1]=1.e+35;
  
/**----------------------------------------------------------------------------
  compute N² at levels : N²=(-g/rho dhro /dz) */
  for (l=1;l<nlevels-1;l++) {
    k=l;
    double r=0.5*(rho[k]+rho[k-1]);
    N2[l]=-g*(rho[k]-rho[k-1])/dZ[l]/r;
    N2[l]=max(N2[l], 1.e-10);
    if(N2[l]<0) {
      printf("compute_Pmode: unstable stratification, abort...\n");
      return(-1);
      }
    }

  N2[0]=N2[1];
  N2[nlevels-1]=N2[nlevels-2];
  
  for(k=0;k<dim*dim;k++) A[k]=0.0;
  for(k=0;k<dim*dim;k++) B[k]=0.0;

/**----------------------------------------------------------------------------
  A matrix, central part, skip the first layer (bottom) and last layer (surface)*/
  for (int k=1;k<nlayers-1;k++) {
    row=k;
    col=row-1;
    A[dim*col+row] = -1./(rho[k]-rho[k-1]);
    col=row;
    A[dim*col+row] =  1./(rho[k]-rho[k-1])+1./(rho[k+1]-rho[k]);
    col=row+1;
    A[dim*col+row] = -1./(rho[k+1]-rho[k]);
    }
    
/**----------------------------------------------------------------------------
  B matrix */
  for (int k=0;k<nlayers;k++) {
    row=k;
    col=row;
    B[dim*col+row] = -g*dz[k];
    }

/**----------------------------------------------------------------------------
  bottom layer*/
  k=0;
  row=k;
  col=row;
  A[dim*col+row] =  1/(rho[k+1]-rho[k]);
  col=row+1;
  A[dim*col+row] = -1/(rho[k+1]-rho[k]);
  
/**----------------------------------------------------------------------------
  top layer */
  k=nlayers-1;
  row=k;
  col=row-1;
  A[dim*col+row] = -1./(rho[k]-rho[k-1]);
  col=row;
  A[dim*col+row] =  1./(rho[k]-rho[k-1])-1./rho[k];
  
  if(debug) status=matrix_print(A,neq,neq);
  if(debug) status=matrix_print(B,neq,neq);
    
/**----------------------------------------------------------------------------
  then solve the eigenvalue problem; imaginary part issue to be discussed */
  char JOBVL='N',JOBVR='V';
  dggev_( &JOBVL, &JOBVR, &neq, A, &dim, B, &dim, ALPHAR, ALPHAI, BETA, VL, &dim, VR, &dim, WORK, &LWORK, &status );
  
  if(verbose==1) STDOUT_BASE_LINE("After dgeev solver: %d  \n",status);
  
  if(status!=0) return(-1);
  
  for (int k=0;k<nlevels-1;k++) {
    for (int j=0;j<nlevels-1;j++) {
      modes[k][j]=0;
      }
    }
  
  nmodes=neq;
  for (int k=0;k<neq;k++) {
    if(fabs(BETA[k]) > 1.e-10) ALPHAR[k]/=BETA[k];
    else {
      ALPHAR[k]=1.e+20;
      nmodes--;
      }
    }
  
/**----------------------------------------------------------------------------
  Enforce decreasing order (not granted by solver) */
  for (int k=0;k<neq;k++) {
/**----------------------------------------------------------------------------
    get the velocity modes from the eigen values vector */
    unsorted[k]=1./sqrt(ALPHAR[k]);
    }
  size_t *order=sort(unsorted,neq);
  
  for (int i=0;i<neq;i++) {
    int k=order[neq-1-i];
    c[i]=1./sqrt(ALPHAR[k]);
    if(ALPHAR[k]==1.e+20) {
      for (int j=0;j<neq+1;j++) {
        modes[i][j]=0;
        }
      continue;
      }
    if(isfinite(c[i])==0) {
      STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
      }
    if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
/**----------------------------------------------------------------------------
    eigenvector from second level to top level */
    for (int j=0;j<neq;j++) {
      if(isfinite(VR[dim*k+j])==0) {
        if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
        }
      modes[i][j]=VR[dim*k+j];
      }
    }
    
  for (int i=0;i<nlevels;i++) {
    N[i]=sqrt(N2[i]);
    }
    
  delete[] dz;
  delete[] dZ;
  delete[] Z;
  delete[] N2;
  
  delete[] order;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_Pmode_v2(int nlevels,double *z, double *rho, double g, double *N, double *c, double **modes, int & nmodes, bool rigid_lid, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  int k,l,row,col;
  int nlayers=nlevels-1;
  double *dZ,*dz,*Z,*N2;
  int dim=nlevels-1,neq=dim;
  
  double A[dim*dim],B[dim*dim];
  int LWORK  =1 + 8*dim;
  double WORK[LWORK],ALPHAR[dim],ALPHAI[dim],BETA[dim],VR[dim*dim],VL[dim*dim];
  double unsorted[dim];
  
  bool debug=false;
  
  dz=new double[nlayers];
  Z =new double[nlayers];
  dZ=new double[nlevels];
  N2=new double[nlevels];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int i=0;i<nlayers;i++) {
    dz[i] = z[i+1]-z[i];
    Z[i]  = 0.5*(z[i+1]+z[i]);
    }

/**----------------------------------------------------------------------------
  Z is the layer depth; dZ is "level" thickness*/
  for (int i=1;i<nlevels-1;i++) {
    dZ[i]=Z[i]-Z[i-1];
    }
  dZ[0]=1.e+35;
  dZ[nlevels-1]=1.e+35;
  
/**----------------------------------------------------------------------------
  compute N² at levels : N²=(-g/rho dhro /dz) */
  for (l=1;l<nlevels-1;l++) {
    k=l;
    double r=0.5*(rho[k]+rho[k-1]);
    N2[l]=-g*(rho[k]-rho[k-1])/dZ[l]/r;
    N2[l]=max(N2[l], 1.e-10);
    if(N2[l]<0) {
      printf("compute_Pmode: unstable stratification, abort...\n");
      return(-1);
      }
    }

  N2[0]=N2[1];
  N2[nlevels-1]=N2[nlevels-2];
  
  for(k=0;k<dim*dim;k++) A[k]=0.0;
  for(k=0;k<dim*dim;k++) B[k]=0.0;

/**----------------------------------------------------------------------------
  A matrix, central part, skip the first layer (bottom) and last layer (surface)*/
  for (int k=1;k<nlayers-1;k++) {
    row=k;
    col=row-1;
    A[dim*col+row] = -1./(rho[k]-rho[k-1]);
    col=row;
    A[dim*col+row] =  1./(rho[k]-rho[k-1])+1./(rho[k+1]-rho[k]);
    col=row+1;
    A[dim*col+row] = -1./(rho[k+1]-rho[k]);
    }
    
/**----------------------------------------------------------------------------
  B matrix */
  for (int k=0;k<nlayers;k++) {
    row=k;
    col=row;
    B[dim*col+row] = -g*dz[k];
    }

/**----------------------------------------------------------------------------
  bottom layer*/
  k=0;
  row=k;
  col=row;
  A[dim*col+row] =  1/(rho[k+1]-rho[k]);
  col=row+1;
  A[dim*col+row] = -1/(rho[k+1]-rho[k]);
  
/**----------------------------------------------------------------------------
  top layer */
  k=nlayers-1;
  row=k;
  col=row-1;
  A[dim*col+row] = -1./(rho[k]-rho[k-1]);
  col=row;
  A[dim*col+row] =  1./(rho[k]-rho[k-1])-1./rho[k];
  
  if(debug) status=matrix_print(A,neq,neq);
  if(debug) status=matrix_print(B,neq,neq);
    
/**----------------------------------------------------------------------------
  then solve the eigenvalue problem; imaginary part issue to be discussed */
  char JOBVL='N',JOBVR='V';
  dggev_( &JOBVL, &JOBVR, &neq, A, &dim, B, &dim, ALPHAR, ALPHAI, BETA, VL, &dim, VR, &dim, WORK, &LWORK, &status );
  
  if(verbose==1) STDOUT_BASE_LINE("After dgeev solver: %d  \n",status);
  
  if(status!=0) return(-1);
  
  for (int k=0;k<nlevels-1;k++) {
    for (int j=0;j<nlevels-1;j++) {
      modes[k][j]=0;
      }
    }
  
  nmodes=neq;
  for (int k=0;k<neq;k++) {
    if(fabs(BETA[k]) > 1.e-10) ALPHAR[k]/=BETA[k];
    else {
      ALPHAR[k]=1.e+20;
      nmodes--;
      }
    }
  
/**----------------------------------------------------------------------------
  Enforce decreasing order (not granted by solver) */
  for (int k=0;k<neq;k++) {
/**----------------------------------------------------------------------------
    get the velocity modes from the eigen values vector */
    unsorted[k]=1./sqrt(ALPHAR[k]);
    }
  size_t *order=sort(unsorted,neq);
  
  for (int i=0;i<neq;i++) {
    int k=order[neq-1-i];
    c[i]=1./sqrt(ALPHAR[k]);
    if(ALPHAR[k]==1.e+20) {
      for (int j=0;j<neq+1;j++) {
        modes[i][j]=0;
        }
      continue;
      }
    if(isfinite(c[i])==0) {
      STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
      }
    if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
/**----------------------------------------------------------------------------
    eigenvector from second level to top level */
    for (int j=0;j<neq;j++) {
      if(isfinite(VR[dim*k+j])==0) {
        if(verbose==1) STDOUT_BASE_LINE("celerity %dth mode : %lf \n",i,c[i]);
        }
      modes[i][j]=VR[dim*k+j];
      }
    }
    
  for (int i=0;i<nlevels;i++) {
    N[i]=sqrt(N2[i]);
    }
    
  delete[] dz;
  delete[] dZ;
  delete[] Z;
  delete[] N2;
  
  delete[] order;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dgeev_interface(double* A, int dim, double** &eigenvectors, double* &eigenvalues, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;

  int neq=dim;
  
  double unsorted[dim];
  int LWORK  =1 + 4*dim;
  double WORK[LWORK],WR[dim],WI[dim],VR[dim*dim],VL[dim*dim];
  

/**----------------------------------------------------------------------------
  then solve the eigenvalue problem; imaginary part issue to be discussed */
  char JOBVL='N',JOBVR='V';
  dgeev_( &JOBVL, &JOBVR, &neq, A, &dim, WR, WI, VL, &dim, VR, &dim, WORK, &LWORK, &status );
  
  if(verbose==1) STDOUT_BASE_LINE("After dgeev solver: %d  \n",status);
  
//   if(status!=0) return(-1);
  if(status!=0) TRAP_ERR_EXIT(-1, "eigenvector computation failed\n");
  
  if(eigenvectors==0) {
    eigenvectors=new double *[dim];
    for (int k=0;k<dim;k++) eigenvectors[k]=new double [dim];
    }
  
  if(eigenvalues==0) {
    eigenvalues=new double [dim];
    }
  
  for (int k=0;k<dim;k++) {
    eigenvalues[k]=0;
    for (int j=0;j<dim;j++) {
      eigenvectors[k][j]=0;
      }
    }
  
/**----------------------------------------------------------------------------
  Enforce decreasing order (not granted by solver) */
  for (int k=0;k<neq;k++) {
/**----------------------------------------------------------------------------
    get the velocity modes from the eigen values vector */
    unsorted[k]=WR[k];
    }
  size_t *order=sort(unsorted,neq);
  
  for (int i=0;i<neq;i++) {
    int k=order[neq-1-i];
    eigenvalues[i]=WR[k];
/**----------------------------------------------------------------------------
    eigenvector from second level to top level */
    for (int j=0;j<neq;j++) {
      if(isfinite(VR[dim*k+j])==0) {
        printf("%dth mode : %lf \n",i,eigenvalues[i]);
        }
      eigenvectors[i][j]=VR[dim*k+j];
      }
    }
      
  return(0);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Wmodes_ProjectionMatrix_v1(int nlevels, double *z, double **modes, int nmodes, double *rho, double *N, double *Pmatrix, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------------
  
  build P projection matrix so that : [a] = P x [u] where [a] are the cofficients
  of u on the modal basis
  
  only for full modal basis, normalization dependant
------------------------------------------------------------------------------*/
{
  int i, k, status;
  int nlayers=nlevels-1;
  
  double *dz=new double[nlayers];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (k=0;k<nlayers;k++) {
    dz[k] = z[k+1]-z[k];
    }
  
  for(i=0;i<nlevels;i++) {
    for (k=0;k<nlevels;k++) {
      Pmatrix[i*nlevels+k]=0.0;
      }
    }
  
//   for(i=0;i<nmodes;i++) {
//     double p[2], q[2], r[2], s[2];
//     status=vector_checkNaN(modes[i], nlevels);
//     for (int k=0;k<nlevels-1;k++) {
//       p[0]=N[k];
//       p[1]=N[k+1];
//       q[0]=1.0; q[1]=0.0;
//       Pmatrix[i*nlevels+k]   +=dz[k]*sqrt(rho[k])*fe_integraleLGP1xLGP1xLGP1_1D(p, &modes[i][k], q);
//       q[0]=0.0; q[1]=1.0;
//       Pmatrix[i*nlevels+k+1] +=dz[k]*sqrt(rho[k])*fe_integraleLGP1xLGP1xLGP1_1D(p, &modes[i][k], q);
//       }
//     }
  for(i=0;i<nmodes;i++) {
    status=vector_checkNaN(modes[i], nlevels);
    for (k=1;k<nlevels-1;k++) {
      Pmatrix[k*nlevels+i]=sqrt(rho[k])*N[k]*(dz[k-1]+dz[k])*modes[i][k]/2.0;
      }
    k=0;
    Pmatrix[k*nlevels+i]=sqrt(rho[k])*N[k]*dz[k]*modes[i][k]/2.0;
    k=nlevels-1;
    Pmatrix[k*nlevels+i]=sqrt(rho[k])*N[k]*dz[k-1]*modes[i][k]/2.0;
    }
  
  delete[] dz;
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Wmodes_ProjectionMatrix_v2(int nlevels, double **modes, int nmodes, double *Pmatrix, bool debug, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------------
  
  build P projection matrix so that : [a] = P x [u] where [a] are the cofficients
  of u on the modal basis
  
  find [a] such that it minimizes |M[a] - [u]|² where M is the modes matrix (column)
  
  Least square solution gives:  P=inv(tMM)tM
  
------------------------------------------------------------------------------*/
{
  int i, k, status;
  
  for(k = 0; k < nlevels*nlevels; k++){
    Pmatrix[k] = 0.0;
    }
  
  if(nmodes==0) return(0);
  
  double *M=new double[nmodes*nlevels], *leastsquare=0;

  for(k = 0; k < nmodes*nlevels; k++){
    M[k] = 0.0;
    }
  
  for(i=0;i<nmodes;i++) {
    for (k=0;k<nlevels;k++) {
      M[i*nlevels+k]=modes[i][k];
      }
    }
  
/**-----------------------------------------------------------------------------
  returns nmodes x nlevels lest-square matrix */
  status=matrix_least_square(M, nlevels, nmodes, leastsquare, verbose);
  if(status!=0) TRAP_ERR_EXIT(-1, "least-square matrix construction failed\n");
    
/**-----------------------------------------------------------------------------
  updates the nmodes firts lines of coefficient Pmatrix with lest-square matrix */
  for(i=0;i<nmodes;i++) {
    for (k=0;k<nlevels;k++) {
      Pmatrix[k*nlevels+i]=leastsquare[k*nmodes+i];
      }
    }
  
  if(verbose==1) status=matrix_print(Pmatrix, nlevels, nlevels);
  
/**-----------------------------------------------------------------------------
  verification, to be used in debug mode only */
  if(debug) {
    double *product=0;
    status=matrix_product(leastsquare,M, nmodes,nlevels,nlevels,nmodes,product,0);
    status=matrix_print(product, nmodes, nmodes);
    delete[] product;
    }
  
  delete[] leastsquare;
  delete[] M;
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Umodes_decomposition_v1(int nlevels, double *z, double **modes, int nmodes, double *rho, complex<double> *u, complex<double> *decomposition, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**------------------------------------------------------------------------------
  
  decompose velocities and level into modal contributions
  
  suitable only for full mode basis, basis orthonormality required (use of
  scalar product)
  
------------------------------------------------------------------------------*/
{
  int i;
  int nlayers=nlevels-1;
  complex<double> scalar_product;
  
  double *dz=new double[nlayers];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int k=0;k<nlayers;k++) {
    dz[k] = z[k+1]-z[k];
    }
  
//   verbose=1;
  
  for(i=0;i<nmodes;i++) {
    scalar_product=0;
    for (int k=0;k<nlevels-1;k++) {
      scalar_product+=dz[k]*modes[i][k]*u[k]*sqrt(rho[k]);
      }
    if(verbose==1) printf("U-modes %d, scalar product=%lf\n",i,abs(scalar_product));
    decomposition[i]=scalar_product;
    }
  
  delete[] dz;
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Umodes_decomposition_v2(int nlevels, double **modes, int nmodes, double *u, double *decomposition, bool debug, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------------
  
  decompose velocities and level into modal contributions
  
  suitable for incomplete mode basis, basis orthogonality not required
  
------------------------------------------------------------------------------*/
{
  int status, undetermined, unit=1, ntruemodes;
  double umodal[nlevels], chk[nlevels];
  
  double *Pmatrix=new double[nlevels*nlevels];
  
  for(int n=0;n<nlevels*nlevels;n++) Pmatrix[n]=0.0;
  for(int m=0;m<nmodes;m++) decomposition[m]=0.0;
  for(int m=0;m<nmodes;m++) chk[m]=0.0;

/**-----------------------------------------------------------------------------
  check for non-existing modes */
  undetermined=0;
  for(int m=0;m<nmodes;m++) {
    if (isnan(modes[m][0])==1) {
      undetermined++;
      continue;
      }
    double r=ddot_(&nlevels,modes[m],&unit,modes[m],&unit);
    if (r==0) {
      undetermined++;
      continue;
      }
    }
  
  ntruemodes=nmodes-undetermined;
  if(ntruemodes==0) {
    delete[] Pmatrix;
    return(0);
    }
//   bool debug=false;
  status=Wmodes_ProjectionMatrix_v2(nlevels, modes, ntruemodes, Pmatrix, debug, 0);
  
//   verbose=1;
  if(verbose==1) status=matrix_print(Pmatrix, nlevels, nlevels);
  
  status=matrix_operation(Pmatrix, u, nlevels, nlevels, decomposition, verbose);
  for(int m=0;m<nmodes;m++) {
    double d=decomposition[m];
    if(!isnormal(d) and d!=0) {
      printf("%s anomaly\n",__func__);
      decomposition[m]=0;
      }
    }

/**-----------------------------------------------------------------------------
  verification */
//   for(int m=0;m<ntruemodes;m++) {
//     for(int k=0;k<nlevels;k++) {
//       chk[m]+=Pmatrix[k*nlevels+m]*u[k];
//       }
//     printf("%d %7.4lf %7.4lf %7.4lf %7.4lf\n",m,real(decomposition[m]), imag(decomposition[m]),real(chk[m]), imag(chk[m]));
//     }
  
/**-----------------------------------------------------------------------------
  verification : should find unity decomposition (i.e. umodal[k]=u[k] when u=mode) */
//   debug=true;

  if(debug) {
    double *buffer=new double[nlevels];
    for(int m=0;m<ntruemodes;m++) {
      for(int k=0;k<nlevels;k++) {
        buffer[k]=modes[m][k];
        }
      status=matrix_operation(Pmatrix, buffer, nlevels, nlevels, decomposition, verbose);
      printf("mode %d\n",m);
      for(int k=0;k<nlevels;k++) {
        umodal[k]=0;
        for(int m=0;m<ntruemodes;m++) {
          umodal[k]+=decomposition[m]*modes[m][k];
          }
        printf("%d %7.4lf %7.4lf\n",k,buffer[k], umodal[k]);
        }
      }
    delete[] buffer;
    }
  

  delete[] Pmatrix;
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Umodes_decomposition_v2(int nlevels, double **modes, int nmodes, complex<double> *u, complex<double> *decomposition, bool debug, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------------
  
  decompose velocities and level into modal contributions
  
  suitable for incomplete mode basis, basis orthogonality not required
  
------------------------------------------------------------------------------*/
{
  int status, undetermined, unit=1, ntruemodes;
  complex<double> scalar_product, *umodal=0, *chk=0;
  
  if(nlevels==0) {
    TRAP_ERR_EXIT(-1,"exiting\n");
    }
  umodal=new complex<double>[nlevels];
  chk=new complex<double>[nlevels];
  
  double *Pmatrix=new double[nlevels*nlevels];
  
  for(int n=0;n<nlevels*nlevels;n++) Pmatrix[n]=0.0;
  for(int m=0;m<nmodes;m++) decomposition[m]=0.0;
  for(int m=0;m<nmodes;m++) chk[m]=0.0;

/**-----------------------------------------------------------------------------
  check for non-existing modes */
  undetermined=0;
  for(int m=0;m<nmodes;m++) {
    if (isnan(modes[m][0])==1) {
      undetermined++;
      continue;
      }
    double r=ddot_(&nlevels,modes[m],&unit,modes[m],&unit);
    if (r==0) {
      undetermined++;
      continue;
      }
    }
  
  ntruemodes=nmodes-undetermined;
  if(ntruemodes==0) {
    delete[] Pmatrix;
    return(0);
    }
//   bool debug=false;
  status=Wmodes_ProjectionMatrix_v2(nlevels, modes, ntruemodes, Pmatrix, debug, 0);
  
//   verbose=1;
  if(verbose==1) status=matrix_print(Pmatrix, nlevels, nlevels);
  
  status=matrix_operation(Pmatrix, u, nlevels, nlevels, decomposition, verbose);
  for(int m=0;m<nmodes;m++) {
    double d=real(decomposition[m]);
    if(!isnormal(d) and d!=0) {
      printf("%s anomaly\n",__func__);
      }
    d=imag(decomposition[m]);
    if(!isnormal(d) and d!=0) {
      printf("%s anomaly\n",__func__);
      }
    }

/**-----------------------------------------------------------------------------
  verification */
//   for(int m=0;m<ntruemodes;m++) {
//     for(int k=0;k<nlevels;k++) {
//       chk[m]+=Pmatrix[k*nlevels+m]*u[k];
//       }
//     printf("%d %7.4lf %7.4lf %7.4lf %7.4lf\n",m,real(decomposition[m]), imag(decomposition[m]),real(chk[m]), imag(chk[m]));
//     }
  
/**-----------------------------------------------------------------------------
  verification : should find unity decomposition (i.e. umodal[k]=u[k] when u=mode) */
//   debug=true;;

  if(debug) {
    complex<double> *buffer=new complex<double>[nlevels];
    for(int m=0;m<ntruemodes;m++) {
      for(int k=0;k<nlevels;k++) {
        buffer[k]=modes[m][k];
        }
      status=matrix_operation(Pmatrix, buffer, nlevels, nlevels, decomposition, verbose);
      printf("mode %d\n",m);
      for(int k=0;k<nlevels;k++) {
        umodal[k]=0;
        for(int m=0;m<ntruemodes;m++) {
          umodal[k]+=decomposition[m]*modes[m][k];
          }
        printf("%d %7.4lf %7.4lf %7.4lf %7.4lf\n",k,real(buffer[k]), real(umodal[k]), imag(buffer[k]), imag(umodal[k]));
        }
      }
    delete[] buffer;
    }
  

  delete[] umodal;
  delete[] Pmatrix;
  delete[] chk;

  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Umodes_decomposition_v2(int nlevels, double **modes, int nmodes, complex<float> *u, complex<float> *decomposition, bool debug, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------------
  
  decompose velocities and level into modal contributions
  
  suitable for incomplete mode basis, basis orthogonality not required
  
------------------------------------------------------------------------------*/
{
  int status;
  complex<float> *uu=0, *dd=0;
  
  recast(u,uu,nlevels);
  recast(decomposition,dd,nlevels);
  
  status=Umodes_decomposition_v2(nlevels, modes, nmodes, uu, dd, debug, verbose);
  
  recast(uu,u,nlevels);
  recast(dd,decomposition,nlevels);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Wmodes_decomposition(int nlevels, double *z, double **modes, int nmodes, double *rho, double *N, complex<double> *u, complex<double> *decomposition, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**-----------------------------------------------------------------------------
  only for full modal basis, orthogonality required, normalization dependant
------------------------------------------------------------------------------*/
{
  int i;
  int nlayers=nlevels-1;
  complex<double> scalar_product;
  
  double *dz=new double[nlayers];

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int k=0;k<nlayers;k++) {
    dz[k] = z[k+1]-z[k];
    }
  
//   verbose=1;
       
  for(i=0;i<nmodes;i++) {
    scalar_product=0;
    double p[2];
    for (int k=0;k<nlevels-1;k++) {
      p[0]=N[k];
      p[1]=N[k+1];
      scalar_product+=dz[k]*sqrt(rho[k])*fe_integraleLGP1xLGP1xLGP1_1D(p, &modes[i][k], &u[k]);
      }
    if(verbose==1) printf("W-modes %d, scalar product=%lf\n",i,abs(scalar_product));
    decomposition[i]=scalar_product;
    }
    
  delete[] dz;
   
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Wmodes_decomposition(int nlevels, double *z, double **modes, int nmodes, double *rho, complex<double> *u, complex<double> *decomposition, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, status;
  int nlayers=nlevels-1;
  complex<double> scalar_product;
  double g=9.81;
  
  double *dz=new double[nlayers], *Z=new double[nlayers];
  double *N=new double[nlevels],*dZ=new double[nlevels], N2;

/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int i=0;i<nlayers;i++) {
    dz[i] = z[i+1]-z[i];
    Z[i]  = 0.5*(z[i+1]+z[i]);
    }

/**----------------------------------------------------------------------------
  Z is the layer depth; dZ is "level" thickness*/
  for (int i=1;i<nlevels-1;i++) {
    dZ[i]=Z[i]-Z[i-1];
    }
  dZ[0]=1.e+35;
  dZ[nlevels-1]=1.e+35;
  
/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int k=0;k<nlayers;k++) {
    dz[k] = z[k+1]-z[k];
    }

/**----------------------------------------------------------------------------
  compute N² at levels : N²=(-g/rho dhro /dz) */
  for (l=1;l<nlevels-1;l++) {
    k=l;
    double r=0.5*(rho[k]+rho[k-1]);
    N2=-g*(rho[k]-rho[k-1])/dZ[l]/r;
    N2=max(N2, 1e-12);
    if(N2<0) {
      printf("Wmodes_decomposition: unstable stratification, abort...\n");
      return(-1);
      }
    N[l]=sqrt(N2);
    }

  N[0]=N[1];
  N[nlevels-1]=N[nlevels-2];
  
  status=Wmodes_decomposition(nlevels, z, modes, nmodes, rho, N, u, decomposition, verbose);
  
  delete[] dz;
  delete[] dZ;
  delete[] N;
  delete[] Z;
   
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Wmodes_unresolved(int nlevels, double *z, double **modes, int nmodes, double *rho, double *N, double *coefficients, double *inverse, bool debug, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, m, n, status;
  int nlayers=nlevels-1;
  double **eigenvectors=0, *eigenvalues=0;
  
  double *dz=new double[nlayers];
  double *r=new double[nlevels];
  
  double *A=new double[nlevels*nlevels];
  
  
  for (k=1;k<nlevels-1;k++) {
    r[k]=0.5*(rho[k]+rho[k-1]);
    }
  r[0]=rho[0];
  r[nlevels-1]=rho[nlayers-1];
  
  
/**----------------------------------------------------------------------------
  z is the level depth; dz is layer thickness*/
  for (int k=0;k<nlayers;k++) {
    dz[k] = z[k+1]-z[k];
    }

/*------------------------------------------------------------------------------
  Diagnostic block : verification of mode's orthogonality */
  for(n=0;n<nlevels*nlevels;n++) A[n]=0.0;
  for(l=0;l<nmodes;l++) {
    for (int k=0;k<nmodes;k++) {
      n=l*nmodes+k; /* fortran ordering  : row k, column l */
      for (int i=0;i<nlevels;i++) {
        m=l*nlevels+i; /* fortran ordering  : row i, column l */
//         M[n]+=coefficients[m]*modes[l][i]*sqrt(r[i])*N[i];
        A[n]+=coefficients[m]*modes[l][i];
        }
      }
    }
    
  if(verbose==1) status=matrix_print(A, nlevels, nlevels);
  
/*------------------------------------------------------------------------------
  compute eigenvalue matrix for additional conditions using modal solution constraint */
  for(n=0;n<nlevels*nlevels;n++) A[n]=0.0;
  for(l=0;l<nlevels;l++) {
    for (int k=0;k<nlevels;k++) {
      n=l*nlevels+k; /* fortran ordering  : row k, column l */
      for (int i=0;i<nmodes;i++) {
/*------------------------------------------------------------------------------
        coefficient[l*mesh.nlevels+k] = coefficient[level k for mode l] */
        m=l*nlevels+i; /* standard ordering  : row i, column l */
        A[n]-=modes[i][k]*coefficients[m];
        }
      }
    }
    
  for(l=0;l<nlevels;l++) {
    n=l*nlevels+l; /* row l, column l */
    A[n]+=1.0;
    }
    
  status=dgeev_interface(A, nlevels, eigenvectors, eigenvalues, verbose);
  
  if(verbose==1) {
    for(l=0;l<nlevels;l++) {
      printf("%3d : %12.6lf -> ", l, eigenvalues[l]);
      for (int k=0;k<nlevels;k++) {
        printf("%12.6lf ",eigenvectors[l][k]);
        }
      printf("\n");
      }
    }
  
//   debug=true;
  if(debug) {
    vector<int> nzero;
    vector<int> nunit;
    for(l=0;l<nlevels;l++) {
      if(fabs(eigenvalues[l])     <1.e-12) nzero.push_back(l);
      if(fabs(1.0-eigenvalues[l]) <1.e-12) nunit.push_back(l);
      }
    printf("nlevels=%d, nzero=%d, nunit=%d\n", nlevels, nzero.size(), nunit.size());
    for(k=0;k<nunit.size();k++) printf("%d ", nunit[k]);
    printf("\n");
    for(k=0;k<nzero.size();k++) printf("%d ", nzero[k]);
    printf("\n");
    }
    
/*------------------------------------------------------------------------------
  bottom level (and surface level if rigid lid) will have a special shape */
  int top=-1,bottom=-1;
  
  for(l=0;l<nlevels;l++) {
    int nzero=0;
    int nunit=0;
    int pos=-1;
/*------------------------------------------------------------------------------
    count number of zeroes and 1 in the l-th eigenvector*/
    for (int k=0;k<nlevels;k++) {
      if(eigenvectors[l][k]==0.0) nzero++;
      if(eigenvectors[l][k]==1.0) {
        pos=k;
        nunit++;
        }
      }
    if( (nunit==1) && (nzero ==nlevels-1) ) {
      if(pos==0) {
        if(verbose==1) printf("bottom mode at %d\n",l);
        bottom=l;
        }
      if(pos==nlevels-1) {
        if(verbose==1) printf("top mode at %d\n",l);
        top=l;
        }
      }
    }
    
  double *tmp;
  int pos=nlevels-nmodes-1;
 
  if(bottom!=-1) {
    tmp=eigenvectors[pos];
    eigenvectors[pos]=eigenvectors[bottom];
    eigenvectors[bottom]=tmp;
    pos--;
    }
    
  if(top!=-1) {
    tmp=eigenvectors[pos];
    eigenvectors[pos]=eigenvectors[top];
    eigenvectors[top]=tmp;
    pos--;
    }
      
/*------------------------------------------------------------------------------
  eigenvector matrix, to be inverted */
  double *Q=new double[nlevels*nlevels];
//   double *inverse=0;
  
  for(l=0;l<nlevels;l++) {
    for (int k=0;k<nlevels;k++) {
      n=l*nlevels+k; /* row k, column l */
      Q[n]=eigenvectors[l][k];
      }
    }
  status=matrix_inverse(Q, nlevels, inverse, false, verbose);

  if(verbose==1) status=matrix_print(inverse, nlevels, nlevels);
  
  delete[] dz;
  delete[] r;
  delete[] A;
  for(l=0;l<nlevels;l++) {
    delete[] eigenvectors[l];
    }
  delete[] eigenvectors;
  delete[] eigenvalues;
  delete[] Q;
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int debug_vertical_mode()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  testing vertical normal modes computation

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@**/
{
  int i,k,l,status;
  int nlayers=20,nlevels=nlayers+1;
  double *z, *Z, *N, *rho, *c, **modes, g=9.81;
  double depth=1000.,dz=depth/nlayers;
  double r;
  FILE *out;
  
  z     =new double[nlevels];
  Z     =new double[nlevels];
  N     =new double[nlevels];
  rho   =new double[nlayers];
  c     =new double[nlevels];
  modes =new double*[nlevels];
  
/**----------------------------------------------------------------------------
  initialize levels */
  for(l=0;l<nlevels;l++) {
    z[l]=-depth+l*dz;
    modes[l] =new double[nlevels];
    }

/**----------------------------------------------------------------------------
  initialize density for uniform N*/

//  r=water_density(10., 35., -500.); TODO
  r=1025;
  for(k=0;k<nlayers;k++) {
    l=k;
    Z[k]=0.5*(z[l+1]+z[l]);
//    rho[k]=r*exp((Z[k]/Z[0]-1.)*2.e-03*2.e-03/g);
    rho[k]=r*exp(-Z[k]*2.e-03*2.e-03/g);
    }
    
  for(l=0;l<nlevels;l++) {
    N[l]=-1.;
    }
    
/**----------------------------------------------------------------------------
  solve for normal modes */
  status= compute_Wmode_v2(nlevels, z,  rho, g, N, c, modes);
  
  out=fopen("normal-modes","w");
  for(l=0;l<nlevels;l++) {
    fprintf(out,"%lf %lf",z[l],N[l]);
    for(i=0;i<nlevels;i++) fprintf(out," %lf",modes[i][l]);
    fprintf(out,"\n");
    }
  fclose(out);
  
  delete[] z;
  delete[] Z;
  delete[] N;
  delete[] rho;
  delete[] c;
  for(l=0;l<nlevels;l++) delete[] modes[l];
  delete[] modes;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int compute_vertical_mode(mesh_t mesh,int pair)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@**/
{
  int status;
  char filename[1024];
//   int layers=mesh.nlayers;
  discretisation_t z_descriptor, u_descriptor;
  int z_discretisation, u_discretisation;

  status=paire_discretisation_id(pair, &z_discretisation, &u_discretisation);
  z_descriptor=get_descriptor(mesh,z_discretisation);

//   int znodes=z_descriptor.nnodes;

//   double Z[layers], density[layers], *modes[znodes];

/*-----------------------------------------------------------------------------
  get mid-layer depth*/
//   for (int i=0;i<z_descriptor.nnodes;i++){
//      for (int k=0;k<layers;k++){
//       Z[i][k]= (z_descriptor.nodes[n].zlevel[k+1]+z_descriptor.nodes[n].zlevel[k])/2.;
//       density[i][k]=1.;
//       }
//     }
    
  printf("in compute_vertical_mode\n");
//   for (int i=0;i<z_descriptor.nnodes;i++){ TODO
//     modes[i]=new double[layers];
//     for (int k=0;k<layers;k++){
//       Z[k]= (z_descriptor.nodes[i].zlevel[k+1]+z_descriptor.nodes[i].zlevel[k])/2.;
//       density[k]=1.;
//       modes[i][k]=0.0;
//       }
//     printf(" vmode over node: %d  \n",i);
//     if (i==0) status=compute_column_vmode_v1(layers,Z,density,modes[i]);
//     }
  sprintf(filename,"troulala.nc");
//   status=archiving_UGdummy3D(filename, mesh, "vmodes", modes,(int) LGP0,(int) LAYERS, (int) LGP0);
  return(0);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int BruntVassala(grid_t & grid, float *R, float fmask, grid_t & topogrid, float *topo, float topomask, int i, int j, float *N, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,l,m,count=0;
  double x,y;
  double rho[grid.nz],Z[grid.nz],zz[grid.nz],dz[grid.nz],N2[grid.nz+1];
  double g=9.81, mask=fmask;
  int nlevels;
  float h;

  m=j*grid.nx+i;
/**----------------------------------------------------------------------------
  get true depth */
  x=map_grid_x (grid, i, j);
  y=map_grid_y (grid, i, j);
  x=map_recale(topogrid,x);
  status=map_interpolation(topogrid, topo,topomask,x,y,&h);
  
/**----------------------------------------------------------------------------
  build (reversed) density profile */
  for(k=0;k<grid.nz;k++) {
    l=grid.nz-1-k;
    size_t n=k*grid.ny*grid.nx+m;
    if(R[n]!=fmask) {
      rho[l]=R[n];
      }
    else {
      rho[l]=mask;
      }
    if(rho[l] != mask) count++;
    Z[l]=-grid.z[k*grid.ny*grid.nx+m];
    }

/**----------------------------------------------------------------------------
  density given by layers, N by levels */
  nlevels=grid.nz+1;
  for(k=0;k<grid.nz-1;k++) {
    dz[k]=Z[k+1]-Z[k];
    l=k+1;
    zz[l]=0.5*(Z[k+1]-Z[k]);
    if( (rho[k]!=mask) && (rho[k+1]!=mask) ) {
      double r=0.5*(rho[k]+rho[k+1]);
      N2[l]=-g*(rho[k+1]-rho[k])/dz[k]/r;
      if(N2[l]<0.) {
//       printf("unstable stratification %lf, abort...\n",N2[l]);
        N2[l]=0.e+00;
        }
      }
    else {
      N2[l]=mask;
      }
    }

  if(count!=0) {
//     N2[0]=0.;
//     N2[nlevels-1]=0.;
    N2[0]=N2[1];
    N2[nlevels-1]=N2[nlevels-2];
    }
  else {
    N2[0]=mask;
    N2[nlevels-1]=mask;
    }
  
  zz[0]=h;
  zz[nlevels-1]=0.;
  
//   for(l=0;l<nlevels;l++) {
//     size_t n=l*grid.ny*grid.nx+m;
  for(k=0;k<nlevels;k++) {    /// HERE !!!
    l=grid.nz-1-k;
    size_t n=k*grid.ny*grid.nx+m;
    if(N2[l]!=mask) {
      N[n]=sqrt(N2[l]);
      if(N[n]>5.0e-02) {
        printf("%d %f\n",l,N[n]);
        }
      }
    else {
      N[n]=fmask;
      }
    z[k*grid.ny*grid.nx+m]=zz[l];
//    z[l*grid.ny*grid.nx+m]=zz[l];
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int BruntVassala(grid_t & grid, grid_t & w_grid, float *T, float *S, float *R, float fmask, grid_t & topogrid, float *topo, float topomask, int i, int j, float *N, double *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  
  T : potential temperature, ° Celsius (reference=?)
  
  S : salinity, PSU
  
  Definition for temperature, salinity and density:
  
  http://journals.ametsoc.org/doi/abs/10.1175/1520-0485%281997%29027%3C0237:ANDVFT%3E2.0.CO;2
  
  Brunt-Vaisala computation:
  
  ftp://io.usp.br/lado_20130322/papers/houry_etal_1987.pdf
  
----------------------------------------------------------------------------*/
{
  int status;
  int k,l,m,count=0;
  double x,y;
  const int nlayers=grid.nz, nlevels=nlayers+1;
  double rho[nlayers],temperature[nlayers],salinity[nlayers];
  double z_layer[nlayers];
  double dz_layer[nlevels],N2[nlevels],z_level[nlevels];
  double g=9.81, mask=fmask;
  float h;
  int unstable=0;
  
  m=j*grid.nx+i;
  
/**----------------------------------------------------------------------------
  get true depth if needed */
  if(topo!=0) {
    x=map_grid_x (grid, i, j);
    y=map_grid_y (grid, i, j);
    x=map_recale(topogrid,x);
    status=map_interpolation(topogrid, topo,topomask,x,y,&h);
    }
  
  bool ascending;
  double factor;
  
  status=check_vertical_direction(grid,&ascending,&factor);
  
/**----------------------------------------------------------------------------
  keep negative depths */
  if(factor==-1) {
    ascending=!ascending;
    factor=1.0;
    }
/**----------------------------------------------------------------------------
  force negative depths */
  else {
    ascending=!ascending;
    factor=-1.0;
    }
  
  if(ascending) {
/**----------------------------------------------------------------------------
    keep natural (ascending) density profile */
    for(k=0;k<nlayers;k++) {
      l=k;
      size_t n=k*grid.ny*grid.nx+m;
      if(R[n]!=fmask) {
        rho[l]         =R[n];
        if(T!=0) temperature[l] =T[n];
        if(S!=0) salinity[l]    =S[n];
        }
      else {
        rho[l]         =mask;
        temperature[l] =mask;
        salinity[l]    =mask;
        }
      if(rho[l] != mask) count++;
      z_layer[l]=factor*grid.z[k*grid.ny*grid.nx+m];
      }
    }
  else {
/**----------------------------------------------------------------------------
    build reversed (ascending) density profile */
    for(k=0;k<nlayers;k++) {
      l=nlayers-1-k;
      size_t n=k*grid.ny*grid.nx+m;
      if(R[n]!=fmask) {
        rho[l]         =R[n];
        if(T!=0) temperature[l] =T[n];
        if(S!=0) salinity[l]    =S[n];
        }
      else {
        rho[l]         =mask;
        temperature[l] =mask;
        salinity[l]    =mask;
        }
      if(rho[l] != mask) count++;
      z_layer[l]=factor*grid.z[k*grid.ny*grid.nx+m];
      }
    }

/**----------------------------------------------------------------------------
  density given by layers, N by levels */
  for(l=0;l<nlevels;l++) N2[l]=mask;
  
  for(k=0;k<nlayers-1;k++) {
    dz_layer[k+1]=z_layer[k+1]-z_layer[k];
    }
    
  if(w_grid.z==0) {
    for(l=1;l<nlevels-1;l++) {
      k=l-1;
      z_level[l]=0.5*(z_layer[k+1]+z_layer[k]);
      }
    if(topo!=0)
      z_level[0]=min((double) h,z_level[1]-dz_layer[1]);
    else
      z_level[0]=z_level[1]-dz_layer[1];
    z_level[nlevels-1]=0;
    }
  else {
    if(ascending) {
      for(l=0;l<nlevels;l++) {
        size_t n=l*grid.ny*grid.nx+m;
        z_level[l]=factor*w_grid.z[n];
        }
      }
    else {
      for(l=0;l<nlevels;l++) {
        size_t n=l*grid.ny*grid.nx+m;
        z_level[nlevels-l-1]=factor*w_grid.z[n];
        }
      }
    }
    
  if(count>1)
  for(l=1;l<nlevels-1;l++) {
    k=l-1;
    if( (rho[k]!=mask) && (rho[k+1]!=mask) ) {
/**----------------------------------------------------------------------------
      http://servforge.legi.grenoble-inp.fr/projects/CDFTOOLS/browser/trunk/cdfbn2.f90 */
      double r=0.5*(rho[k]+rho[k+1]);
      double chk;
      chk=-g*(rho[k+1]-rho[k])/dz_layer[l]/r;
      if(chk<0.) {
        chk=0.e+00;
        }
/**----------------------------------------------------------------------------
      compute ambient temperature at layer level */
//       double t1=physic_theta(salinity[k],   temperature[k],   10.0, -z_layer[k]);
//       double t2=physic_theta(salinity[k+1], temperature[k+1], 10.0, -z_layer[k+1]);
      double r1,r2, t1, t2, rlocal, drdz_1;
      if(T!=0 and S!=0) {
/**----------------------------------------------------------------------------
        use potential temperature at layer level */
        t1=temperature[k];
        t2=temperature[k+1];

/**----------------------------------------------------------------------------
        compute ambient density at layer level */
        r1=water_density(t1, salinity[k],   -z_layer[k],   1);
        r2=water_density(t2, salinity[k+1], -z_layer[k+1], 1);
        rlocal=0.5*(r1+r2);
      
        drdz_1=(r2-r1)/dz_layer[l];
/**----------------------------------------------------------------------------
        remove compressibility effects */
//       double t1_a=physic_theta(salinity[k],   temperature[k],   10.0, -z_layer[k+1]);
        double t1_a=temperature[k];
        double r1_a=water_density(t1_a, salinity[k], -z_layer[k+1], 1);
        double drdz_2=(r1_a-r1)/dz_layer[l];
/**----------------------------------------------------------------------------
        remove compressibility effects */
//       double t2_a=physic_theta(salinity[k+1],   temperature[k+1],   10.0, -z_layer[k]);
        double t2_a=temperature[k+1];
        double r2_a=water_density(t2_a, salinity[k+1], -z_layer[k], 1);
        double drdz_3=(r2-r2_a)/dz_layer[l];
        N2[l]=-g*(drdz_1-0.5*drdz_2-0.5*drdz_3)/rlocal;
        }
      else {
        r1=rho[k];
        r2=rho[k+1];
        rlocal=0.5*(r1+r2);
      
        drdz_1=(r2-r1)/dz_layer[l];
        
        N2[l]=-g*drdz_1/rlocal;
        }
      
      if(N2[l]<0.) {
/**----------------------------------------------------------------------------
        compute ambient density at central level */
//         r1=water_density(t1, salinity[k],   -z_level[l], 1);
//         r2=water_density(t2, salinity[k+1], -z_level[l], 1);
//         rlocal=0.5*(r1+r2);
//         chk=-g*(r2-r1)/dz_layer[l]/rlocal;
        N2[l]=1.e-04;
        unstable++;
        }
      }
    else {
      N2[l]=mask;
      }
    }

  if(count>1) {
    N2[0]=N2[1];
    N2[nlevels-1]=N2[nlevels-2];
    int nonzero=0;
    for(k=0;k<nlevels;k++) {
      if(N2[k]!=0.0 && N2[k]!=mask) {
        nonzero++;
        }
      }
    if(nonzero==0) {
      STDOUT_BASE_LINE("%s : anomaly at m=%d (unstable column) i=%d j=%d \n",__func__,m,i,j);
      }
    }
  else {
    N2[0]=mask;
    N2[nlevels-1]=mask;
    }
  
  for(k=0;k<nlevels;k++) {
    if(ascending) {
      l=k;
      }
    else {
      l=nlevels-1-k;
      }
    size_t n=k*grid.ny*grid.nx+m;
    if(N2[l]!=mask) {
      N[n]=sqrt(N2[l]);
      if(N[n]>5e-2) {
        //printf("%s : anomaly at m=%d, level=%2d z=%f N=%lf>5e-2\n",__func__,m,l,z_level[l],N[n]);
        }
      }
    else {
      N[n]=fmask;
      }
    z[k*grid.ny*grid.nx+m]=z_level[l];
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SG_Wmodes_v1(const grid_t & grid, float *R, float mask, grid_t & topogrid, float *topo, float topomask, int i, int j, float *Nbar, float **celerity, float **modes, int nmodes, double dz)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  
  made for layer-averaged density field (most of climatologies)
  
 ----------------------------------------------------------------------------*/
{
  int status;
  int k,l,m;
  double x,y;
  double base_rho[grid.nz],base_z[grid.nz],zmin,zmax;
  double *rho,*N=0,*c,**local_modes,g=9.81,*z;
  double w,sum;
  int nlevels, base_nlevels;
  float h;

/**----------------------------------------------------------------------------
  grid point index */
  m=j*grid.nx+i;
  
/**----------------------------------------------------------------------------
  get true depth */
  x=map_grid_x (grid, i, j);
  y=map_grid_y (grid, i, j);
  x=map_recale(topogrid,x);
  status=map_interpolation(topogrid, topo,topomask,x,y,&h);
  
  bool ascending;
  double factor;
  
  status=check_vertical_direction(grid,&ascending,&factor);
  
/**----------------------------------------------------------------------------
  keep negative depths */
  if(factor==-1) {
    ascending=!ascending;
    factor=1.0;
    }
/**----------------------------------------------------------------------------
  force negative depths */
  else {
    ascending=!ascending;
    factor=-1.0;
    }

  if(ascending) {
    for(k=0;k<grid.nz;k++) {
      l=k;
      base_rho[l]=R[k*grid.ny*grid.nx+m];
      base_z[l]=factor*grid.z[k*grid.ny*grid.nx+m];
      }
    }
  else {
/**----------------------------------------------------------------------------
    build (reversed sign and order) density profile (ORCA and others)
    base array are oriented upward, z is set as negative
    first value at ocean bottom (masked if below bathymetry)*/
    for(k=0;k<grid.nz;k++) {
      l=grid.nz-1-k;
      base_rho[l]=R[k*grid.ny*grid.nx+m];
      base_z[l]=factor*grid.z[k*grid.ny*grid.nx+m];
      }
    }
    
  base_nlevels=0;
  for(k=0;k<grid.nz;k++) {
    if(base_rho[k]!=mask) base_nlevels++;
    }

  if(base_nlevels>2) {

/**----------------------------------------------------------------------------
    first unmasked level at grid.nz-base_nlevels */
    zmin=base_z[grid.nz-base_nlevels];
    zmax=base_z[grid.nz-1];
    
    if(h<zmin) {
/**----------------------------------------------------------------------------
      model level is too deep, correct it (to be fixed) */
      zmin=h;
      base_z[grid.nz-base_nlevels]=h;
      }
      
/**----------------------------------------------------------------------------
    local number of valid levels (therefore local number of modes) */
    nlevels=(zmax-zmin)/dz+1;
    dz=(zmax-zmin)/(nlevels-1);
    rho=new double[nlevels];
    z  =new double[nlevels];
    for(k=0;k<nlevels;k++) {
      z[k]=zmin+k*dz;
      }
    for(k=0;k<nlevels-1;k++) {
/**----------------------------------------------------------------------------
      interpolate rho at layer interface */
      status=map_interpolate1D(&(base_rho[grid.nz-base_nlevels]), &(base_z[grid.nz-base_nlevels]),(double) mask, base_nlevels, 0.5*(z[k]+z[k+1]), &rho[k]);
      if(rho[k]==mask) {
        status=map_interpolate1D(base_rho, base_z,(double) mask, grid.nz, 0.5*(z[k]+z[k+1]), &rho[k]);
        }
      }
      
/**----------------------------------------------------------------------------
    allocate arrays */
    c=new double[nlevels];
    for(l=0;l<nlevels;l++) {
      c[l] =1.0;
      }
      
    N=new double[nlevels];
    
    local_modes =new double*[nlevels];
    for(l=0;l<nlevels;l++) {
      local_modes[l] =new double[nlevels];
      }
      
/**----------------------------------------------------------------------------
    set arrays to mask value */
    for(k=0;k<nmodes;k++) celerity[k][m]=mask;
    Nbar[m]=mask;
    
/**----------------------------------------------------------------------------
    solve for normal modes */
    if(nlevels>2) {
      status= compute_Wmode_v2(nlevels, z,  rho, g, N, c, local_modes);
      }

    if(status==0) {
      for(k=0;k<min(nmodes,nlevels-1);k++) celerity[k][m]=c[k];
      for(k=0;k<min(nmodes,nlevels-1);k++) {
        for(i=0;i<grid.nz;i++) {
          l=grid.nz-1-i;
          double dum,factor=1.;
          status=map_interpolate1D(local_modes[k], z,(double) mask, nlevels, base_z[l], &dum);
          if(local_modes[k][nlevels-1]<0) factor=-1;
          if(dum!=mask) modes[k][i*grid.ny*grid.nx+m]=dum*factor;
          }
        }
//       for(k=0;k<min(nmodes,nlevels-1);k++) {
//         for(l=0;l<nlevels;l++) {
//           modes[k][l+grid.nz-nlevels]=local_modes[k][l];
//           }
//         }
/**----------------------------------------------------------------------------
      compute 1st mode-averaged N */
      Nbar[m]=0;
      sum=0;
      for(l=0;l<nlevels;l++) {
        w=local_modes[0][l];
        Nbar[m]+=w*N[l];
        sum+=w;
        }
      Nbar[m]/=sum;
      }
    else {
      status=-1;
      }
    delete[] rho;
    delete[] z;
    delete[] c;
    delete[] N;
    for(l=0;l<nlevels;l++) {
      delete[] local_modes[l];
      }
    delete[] local_modes;
    }
  else {
    status=-1;
    for(k=0;k<nmodes;k++) celerity[k][m]=mask;
    Nbar[m]=mask;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int SG_Wmodes_v2(const grid_t & grid, bool ascending, double factor, float *T, float *S, float *R, float mask, grid_t & topogrid, float *topo, float topomask, int i, int j, float *Nbar, float **celerity, float **modes, int *nmodes, int nkeep,bool rigid_lid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/ 
/**----------------------------------------------------------------------------
  
  made for layer-averaged density field (most of numerical models)
  
  input grid will be assumed to be w-grid
  
 ----------------------------------------------------------------------------*/
{
  int status;
  int k,l,m, nvalues, start, end;
  double x,y;
//   double zmin,zmax;
  double *N=0,*c,**local_modes,g=9.81;
  double w,sum;
  const int nlayers=grid.nz-1, nlevels=nlayers+1;
//   double z_layer[nlayers];
  double *vT=0, *vS=0;
  double vR[nlayers];
  double z_levels[nlevels];
  float h;
  
  if(T!=0 and S!=0) {
    vT=new double[nlayers];
    vS=new double[nlayers];
    }

/**----------------------------------------------------------------------------
  grid point index */
  m=j*grid.nx+i;
  
/**----------------------------------------------------------------------------
  get true depth */
  x=map_grid_x (grid, i, j);
  y=map_grid_y (grid, i, j);
  x=map_recale(topogrid,x);
  status=map_interpolation(topogrid, topo,topomask,x,y,&h);
  
  nvalues=0;
  start=-1;
  end=-1;
  
/**----------------------------------------------------------------------------
  keep negative depths */
  if(factor==-1) {
    ascending=!ascending;
    factor=1.0;
    }
/**----------------------------------------------------------------------------
  force negative depths */
  else {
    ascending=!ascending;
    factor=-1.0;
    }
  
  if(ascending) {
    for(k=0;k<nlayers;k++) {
      if(R[k*grid.ny*grid.nx+m]!=mask) {
        nvalues++;
        if(start==-1) start=k;
        }
      }
    end=nlayers;
    }
  else {
    for(k=0;k<nlayers;k++) {
      if(R[k*grid.ny*grid.nx+m]!=mask) {
        nvalues++;
        if(start==-1) start=k;
        }
      else {
        if(end==-1) end=k;
        }
      }
    start=0;
    if(end==-1) end=nlayers;
    }


  if(nvalues>2) {
    
/**----------------------------------------------------------------------------
  force ascending order (1st level is at bottom) */
  if(ascending) {
    for(k=start;k<end;k++) {
      l=k-start;
      if(T!=0) vT[l]=T[k*grid.ny*grid.nx+m];
      if(S!=0) vS[l]=S[k*grid.ny*grid.nx+m];
      vR[l]=R[k*grid.ny*grid.nx+m];
      }
    for(k=start;k<end+1;k++) {
      l=k-start;
      z_levels[l]=factor*grid.z[k*grid.ny*grid.nx+m];
      }
    }
  else {
/**----------------------------------------------------------------------------
    build (reversed sign and order) density profile (ORCA and others)
    base array are oriented upward, z is set as negative
    first value at ocean bottom (masked if below bathymetry)*/
    for(k=start;k<end;k++) {
      l=end-1-k-start;
      if(T!=0) vT[l]=T[k*grid.ny*grid.nx+m];
      if(S!=0) vS[l]=S[k*grid.ny*grid.nx+m];
      vR[l]=R[k*grid.ny*grid.nx+m];
      }
    for(k=start;k<end+1;k++) {
      l=end-k-start;
      z_levels[l]=factor*grid.z[k*grid.ny*grid.nx+m];
      }
    }


// /**----------------------------------------------------------------------------
//     first unmasked level at grid.nz-base_nlevels */
//     zmin=z_layer[grid.nz-nvalues];
//     zmax=z_layer[grid.nz-1];
    
//     if(h<zmin) {
// /**----------------------------------------------------------------------------
//       model level is too deep, correct it (to be fixed) */
//       zmin=h;
//       z_layer[grid.nz-nvalues]=h;
//       }
          
/**----------------------------------------------------------------------------
    allocate arrays */
    c=new double[nlevels];
    for(l=0;l<nlevels;l++) {
      c[l] =1.0;
      }
      
    N=new double[nlevels];
    local_modes =new double*[nlevels];
    for(l=0;l<nlevels;l++) {
      local_modes[l] =new double[nlevels];
      for(k=0;k<nlevels;k++) local_modes[l][k]=mask;
      }
    
//     z_levels[nlevels-1]=0.5*(z_layer[nlayers-1]-z_layer[nlayers-2])+z_layer[nlayers-1];
//     for(l=nlevels-2; l>=0 ;l--) {
//       z_levels[l]=2*z_layer[l]-z_levels[l+1];
//       }
          
/**----------------------------------------------------------------------------
    set arrays to mask value */
    for(k=0;k<nkeep;k++) celerity[k][m]=mask;
    Nbar[m]=mask;
    
/**----------------------------------------------------------------------------
    solve for normal modes */
//     bool rigid_lid=true;
    int verbose=0;
    status=compute_Wmode_v6(nvalues+1, z_levels, vT, vS, vR, g, N, c, local_modes, nmodes[m], rigid_lid, verbose);

    if(status==0) {
      if(nmodes[m]==0) {
        printf("compute_Wmode_v6: anomaly at m=%d (i=%d j=%d), nmodes=%d\n",m,i,j,nmodes[m]);
        }
      for(int n=0;n<min(nmodes[m], nkeep);n++) {
        celerity[n][m]=c[n];
        double factor=1.;
        if(local_modes[n][nvalues-1]<0) factor=-1;
        for(k=start;k<end+1;k++) {
/**----------------------------------------------------------------------------
          restore original order */
          if(ascending) {
            l=k-start;
            }
          else {
            l=end-k-start;
            }
          modes[n][k*grid.ny*grid.nx+m]=local_modes[n][l]*factor;
          }
        }
/**----------------------------------------------------------------------------
      compute 1st mode-averaged N */
      Nbar[m]=0;
      sum=0;
      for(l=0;l<nmodes[m];l++) {
        w=modes[0][l];
        Nbar[m]+=w*N[l];
        sum+=w;
        }
      Nbar[m]/=sum;
      }
    else {
      status=-2;
      }
    delete[] c;
    delete[] N;
    }
  else {
    status=-1;
    for(k=0;k<nkeep;k++) celerity[k][m]=mask;
    Nbar[m]=mask;
    }
    
  if(vT!=0) delete[] vT;
  if(vS!=0) delete[] vS;
 
  
  return(status);
}

