/** \file
\brief full matrix operation definition
*/

/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/**
\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
  E-mail: florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

***************************************************************************/

#include "config.h"

#include <stdio.h>


#include "tools-define.h"
#include "tools-structures.h"

#include "matrix.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sgmatrix_symmetry(float *A, int dim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n;
  float error=0.;
  int status=0;

  for (l=0;l<dim;l++)
    for (k=l+1;k<dim;k++) {
    m=dim*k+l;
    n=dim*l+k;
    if(A[n]!=A[m])
/*    if(fabs(A[n]-A[m]) > 1.e-06*fabs(A[n]+A[m])) */
      {
/*       printf("gmatrix_symmetry: %d %d %f %f\n",k,l,A[m],A[n]); */
      error=MAX(error,fabs(A[n]-A[m])/fabs(A[n]+A[m]));
      }
    }
  
  if(error !=0.) {
    printf("sgmatrix_symmetry: %g\n",error);
    status=-1;
    }

  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dgmatrix_symmetry(double *A, int dim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n;
  double error=0.;
  int status=0;

  for (l=0;l<dim;l++)
    for (k=l+1;k<dim;k++) {
    m=dim*k+l;
    n=dim*l+k;
    if(A[n]!=A[m])
/*     if(fabs(A[n]-A[m]) > 1.e-06*fabs(A[n]+A[m]))  */
      {
/*       printf("dgmatrix_symmetry: %d %d %lf %lf %lf\n",k,l,A[m],A[n],fabs(A[n]-A[m])/fabs(A[n]+A[m])); */
      error=MAX(error,fabs(A[n]-A[m])/fabs(A[n]+A[m]));
      }
    }
  
  if(error !=0.) {
    printf("dgmatrix_symmetry: %g\n",error);
    status=-1;
    }

  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmatrix_ubw(float *A, int nrow, int ncol)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n,bw=0,l1,l2;

  for(k=0;k<nrow;k++) {
    l1=0;
    l2=ncol-1;
    for(l=0;l<ncol;l++) {
      if(A[l*nrow+k]!=0.) {
        l1=l;
        break;
        }
      }
    if(l==ncol) continue;
    l2=l1;
    for(l=ncol-1;l>l1;l--) {
      if(A[l*nrow+k]!=0.) {
        l2=l;
        break;
        }
      }
    bw=MAX(bw,l2-l1+1);
    }
  return(bw);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmatrix_lbw(float *A, int nrow, int ncol)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n,bw=0,k1,k2;

  for(l=0;l<ncol;l++) {
    k1=0;
    k2=nrow-1;
    for(k=0;k<nrow;k++) {
      if(A[l*nrow+k]!=0.) {
        k1=k;
        break;
        }
      }
    if(k==nrow) continue;
    k2=k1;
    for(k=nrow-1;k>k1;k--) {
      if(A[l*nrow+k]!=0.) {
        k2=k;
        break;
        }
      }
    bw=MAX(bw,k2-k1+1);
    }
  return(bw);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dgmatrix_ubw(double *A, int nrow, int ncol)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n,bw=0,l1,l2;

  for(k=0;k<nrow;k++) {
    l1=0;
    l2=ncol-1;
    for(l=0;l<ncol;l++) {
      if(A[l*nrow+k]!=0.) {
        l1=l;
        break;
        }
      }
    if(l==ncol) continue;
    l2=l1;
    for(l=ncol-1;l>l1;l--) {
      if(A[l*nrow+k]!=0.) {
        l2=l;
        break;
        }
      }
    bw=MAX(bw,l2-l1+1);
    }
  return(bw);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dgmatrix_lbw(double *A, int nrow, int ncol)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,m,n,bw=0,k1,k2;

  for(l=0;l<ncol;l++) {
    k1=0;
    k2=nrow-1;
    for(k=0;k<nrow;k++) {
      if(A[l*nrow+k]!=0.) {
        k1=k;
        break;
        }
      }
    if(k==nrow) continue;
    k2=k1;
    for(k=nrow-1;k>k1;k--) {
      if(A[l*nrow+k]!=0.) {
        k2=k;
        break;
        }
      }
    bw=MAX(bw,k2-k1+1);
    }
  return(bw);

}

/*----------------------------------------------------------------------------*/
///  A is (nrow,ncol), B is (ncol,nrow) */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> int gmatrix_transpose_template(T *A, T *B, int nrow, int ncol)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l;
  for(l=0;l<ncol;l++) {
    for(k=0;k<nrow;k++) {
      B[k*ncol+l]=A[l*nrow+k]; /* A[k,l], B[l,k] */
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dgmatrix_transpose(double *A, double *B, int nrow, int ncol)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=gmatrix_transpose_template(A,B,nrow,ncol);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gmatrix_transpose(float *A, float *B, int nrow, int ncol)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=gmatrix_transpose_template(A,B,nrow,ncol);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *gmatrix_product(float *A, float *B, int dim1, int dim2, int dim3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  A  is square (dim1,dim2)
  B  is square (dim2,dim3)
  C  is square (dim1,dim3)
 */

{
  float *C=NULL,*b=NULL;
  int k,l,m,n;

  exitIfNull(
    C=(float *) malloc(dim1*dim3*sizeof(float))
    );
  exitIfNull(
    b=(float *) malloc(dim2*sizeof(float))
    );

  for(k=0;k<dim1;k++) {
    for(m=0;m<dim2;m++) b[m]=A[m*dim1+k];
    for(l=0;l<dim3;l++) {
#if ATLAS == 1
       C[l*dim1+k]=cblas_dsdot(dim2, b, 1, &(B[l*dim2]), 1);
#else
//       genere une error
#endif
      }
    }

  free(b);
  return(C);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *dmatrix_product(double *A, double *B, int dim1, int dim2, int dim3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  A  is  (dim1,dim2)
  B  is  (dim2,dim3)
  C  is  (dim1,dim3)
 */

{
  double *C=NULL,*b=NULL;
  int k,l,m,n;

  exitIfNull(
    C=(double *) malloc(dim1*dim3*sizeof(double))
    );
  exitIfNull(
    b=(double *) malloc(dim2*sizeof(double))
    );

  for(k=0;k<dim1;k++) {
    for(m=0;m<dim2;m++) b[m]=A[m*dim1+k];
    for(l=0;l<dim3;l++) {
/*
      C[l*dim1+k]=0.;
      for(m=0;m<dim2;m++) {
        C[l*dim1+k]+=A[m*dim1+k]*B[l*dim2+m];
        }
*/

#if ATLAS == 1
C[l*dim1+k] =cblas_ddot(dim2, b, 1,(double *) &(B[l*dim2]), 1) ;
#else
//   genere une error
#endif
      }
    }

  free(b);
  return(C);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *dmatrix_product02(double *A, double *B, int dim1, int dim2, int dim3)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  A  is  (dim1,dim2)
  B  is  (dim2,dim3)
  C  is  (dim1,dim3)
 */

{
  double *C=NULL,*b=NULL;
  int k,l,m,n;
  int *index=NULL,nindex;

  exitIfNull(
    C=(double *) malloc(dim1*dim3*sizeof(double))
    );
  exitIfNull(
    b=(double *) malloc(dim2*sizeof(double))
    );
  exitIfNull(
    index=(int *) malloc(dim2*sizeof(double))
    );

  for(k=0;k<dim1;k++) {
/*-----------------------------------------------------------------------------
    get row k of A */
    nindex=0;
    for(m=0;m<dim2;m++) {
      b[m]=A[m*dim1+k];
      if(b[m] !=0. ) {
        nindex++;
        index[nindex-1]=m;
        }
      }
    if(nindex==0) {
      for(l=0;l<dim3;l++) C[l*dim1+k]=0;
      }

    if(nindex < 100) {
      for(l=0;l<dim3;l++)
/*-----------------------------------------------------------------------------
      scalar product with col l of G */
        {
        C[l*dim1+k]=0;
        for (n=0;n<nindex;n++) {
          m=index[n];
          C[l*dim1+k]+= b[m]*B[l*dim2+m];
          }
        }
      }
    else
      {
      for(l=0;l<dim3;l++)
        #if ATLAS == 1
        C[l*dim1+k]=cblas_ddot(dim2, b, 1,(double *) &(B[l*dim2]), 1);
        #else
//        genere une error
          continue;
        #endif
      }
    }

  free(b);
  free(index);
  return(C);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *bmatrix_product(float *A, float *B, int dim1, int dim2, int dim3, int ldaA, int ldaB)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  A  is square (dim1,dim2) , passed as (ldaA,dim2)
  B  is square (dim2,dim3) , passed as (ldaB,dim3)
  C  is square (dim1,dim3) , passed as (ldaC,dim3)
 */

{
  float *C=NULL;
  int k,l,m,n,r,s;
  int hbwA,hbwB,hbwC;

  hbwA=ldaA/2;
  hbwB=ldaB/2;

  exitIfNull(
    C=(float *) malloc(dim1*dim3*sizeof(float))
    );

  for(l=0;l<dim3;l++) {
    for(k=0;k<dim1;k++) {
      C[l*dim1+k]=0.;
      for(m=0;m<dim2;m++) {
          r=0;
        C[l*dim1+k]+=A[m*dim1+k]*B[l*dim2+m];
        /*       C[l*dim1+k]+=A[m*dim1+k]*B[l*dim2+m];*/
        }
      }
    }

  return(C);

}
