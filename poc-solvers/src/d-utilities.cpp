
/**************************************************************************

  POC-SOLVERS interface, 2006-2012

  Part of the Unstructured Ocean Grid initiative and SYMPHONIE suite

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

***************************************************************************/

#define VERBOSE

#include <string.h>
#include <omp.h>
#include "solverlib.h"
#include "cs.h"

/*=======================================================================*/
/*======================UTILITAIRES MODE PARALLEL =======================*/
/*=======================================================================*/

#ifdef HAVE_MPI

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int d_coo2glob(MPI_Comm world, int ncol, int nnz, int *Ap, int *Ai, double *Ax,
               int *ncolglob, int *nnzglob, int **Apglobin, int **Aiglobin, double **Axglobin) 

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  (Ap, Ai, Ax in COO FORMAT)   
 
  Ax : matrix coefficients
  Ap : rows   vector
  Ai : column vector
  
  still some int/long vulnerability

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{  
  int mpid,sz;
  int nc,nn;
  int *Aptmpg;
  int *Aitmpg;
  double *Axtmpg, *Axglob;
  int bcl;
  int *debtab;
  int *debtab2; 
  int *Apglob;
  int *Aiglob;
  bool debug=false;
    
  MPI_Comm_rank(world, &mpid);
  MPI_Comm_size(world, &sz);
  
/*------------------------------------------------------------------------------
  evaluate global values for number of columns and non-zero values */
  nc=ncol;
  nn=nnz;
  
  MPI_Allreduce(&nc, ncolglob, 1, MPI_INTEGER, MPI_SUM, world);
  MPI_Allreduce(&nn, nnzglob,  1, MPI_INTEGER, MPI_SUM, world);
  MPI_Barrier(world);

  if(debug) printf("%s cpu=%3d : global ndof=%d nnz=%d\n", __func__, mpid, *ncolglob, *nnzglob);

/*------------------------------------------------------------------------------
  allocate global arrays */
  nn = *nnzglob;
  Aptmpg=(int *)malloc(sizeof(int)*nn);
  Apglob=(int *)malloc(sizeof(int)*nn);

  Aitmpg=(int *)malloc(sizeof(int)*nn);
  Aiglob=(int *)malloc(sizeof(int)*nn);

  Axtmpg=(double *)malloc(sizeof(double)*nn);
  Axglob=(double *)malloc(sizeof(double)*nn);

/*------------------------------------------------------------------------------
  initialize global arrays */
  memset(Aptmpg, 0, nn*sizeof(int));
  memset(Aitmpg, 0, nn*sizeof(int));
  memset(Axtmpg, 0.0, nn*sizeof(double));
    
/*------------------------------------------------------------------------------
  each proc fills in his part of global matrices */
  debtab  =(int *)malloc(sizeof(int)*(sz+1));
  debtab2 =(int *)malloc(sizeof(int)*(sz+1));

  memset(debtab,  0, (sz+1)*sizeof(int));
  memset(debtab2, 0, (sz+1)*sizeof(int));

  for (bcl=mpid+1;bcl<(sz+1);bcl++)  {
    debtab2[bcl]=nnz;
    }
  MPI_Allreduce(debtab2, debtab, sz+1, MPI_INTEGER, MPI_SUM, world);
  free(debtab2);
  
  for (bcl=0; bcl < nnz; bcl++) {
    Aptmpg[debtab[mpid]+bcl] = Ap[bcl];
    Aitmpg[debtab[mpid]+bcl] = Ai[bcl];
    Axtmpg[debtab[mpid]+bcl] = Ax[bcl];
    }
  free(debtab);
    
  MPI_Allreduce(Aptmpg, Apglob, nn, MPI_INTEGER, MPI_SUM, world);
  MPI_Allreduce(Aitmpg, Aiglob, nn, MPI_INTEGER, MPI_SUM, world);
  MPI_Allreduce(Axtmpg, Axglob, nn, MPI_DOUBLE , MPI_SUM, world);
  
  free(Aptmpg);
  free(Aitmpg);
  free(Axtmpg);
  
  *Apglobin = Apglob;
  *Aiglobin = Aiglob;
  *Axglobin = Axglob;
    
  return EXIT_SUCCESS;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int d_globcoo2csc(MPI_Comm world, int ncolglob, int nnzglob, int *Apglob,int *Aiglob,double *Axglob,
                    int **Ap,int **Ai,double **Ax)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
 
  global COO vectors   
 
  Axglob : matrix coefficients
  Apglob : rows   vector
  Aiglob : column vector
  
  Ax  : matrix coefficients
  Ap  : rows   vector
  Ai  : column vector
  
  still some int/long vulnerability

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  /* Thanks to Yousef Saad SPARSKIT */
  int mpid,sz;
  int *Aptmpg;
  int *Aitmpg;
  double *Axtmpg;
  int ncol, nnz, k, i, j, k0, pos, bcl, bcl2;
  double x;
  int debug=0;
     
  MPI_Comm_rank(world, &mpid);
  MPI_Comm_size(world, &sz);
  
  ncol=ncolglob;
  nnz=nnzglob;
  
  Aptmpg=(int *)malloc(sizeof(int   )*(ncol+1));
  Aitmpg=(int *)malloc(sizeof(int   )*nnz);
  Axtmpg=(double *)malloc(sizeof(double)*nnz);
    
//   for(bcl=0; bcl<ncol+1;bcl++) Aptmpg[bcl]=0;
  memset(Aptmpg, 0, (ncol+1)*sizeof(int));
  MPI_Barrier(world);

/*------------------------------------------------------------------------------
  Determine the row lengths, use Aptmpg as cardinal counter */
  for(k = 0; k < nnz; k++) {
/*------------------------------------------------------------------------------
    Apglob[k] is column index */
    Aptmpg[Apglob[k]] ++;
    }
  MPI_Barrier(world);

/*------------------------------------------------------------------------------
  The starting position of each row, set Aptmpg as pointer */
  k = 0;
  for (j=0; j<ncol+1; j++) {
    k0 = Aptmpg[j];
    Aptmpg[j] = k;
    k = k+k0;
    }
   
  MPI_Barrier(world);

/*------------------------------------------------------------------------------
  Go through the structure once more.  Fill in output matrix */
  for (k=0;k<nnz;k++) {
    i =  Apglob[k]; // Fortran Numbering (start with 1)
    j =  Aiglob[k]; // Fortran Numbering
    x =  Axglob[k];
    pos = Aptmpg[i];
    Axtmpg[pos] = x;
    Aitmpg[pos] = j; // Fortran Numbering
    Aptmpg[i] = pos+1;
    }
   
  MPI_Barrier(world);
  
/*------------------------------------------------------------------------------
  at that point Aptmpg[i] pointer to incidence of i+1 */
//   Shift back IAO.
//   for (j=ncol; j > 0; j--) {
//     Aptmpg[j] = Aptmpg[j-1];
//     }
/*------------------------------------------------------------------------------
  Fortran trick: i ->i+1 (already done), pos->pos+1 (to do) */
  for (i=0; i<ncol+1;i++) Aptmpg[i]++; // Fortran Numbering
              
  *Ap=Aptmpg;
  *Ai=Aitmpg;
  *Ax=Axtmpg;

  return EXIT_SUCCESS;
}

#endif