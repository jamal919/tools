
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "functions.h"
#include "tools-structures.h"

typedef struct {
   int nx,ny;
   int ncmax,nrmax;
   int *nc,*nr;
   int *row,*col;
   double *matrix;
   } sparse_t ;

#if 0
//unused
/*----------------------------------------------------------------------------*/

int Gsolve(double *A,int neq,int *pivot, double *b)

/*----------------------------------------------------------------------------*/
{
  int status,nrhs=1;

  status=poc_getrf(neq, A, pivot);

  status=poc_getrs(neq,nrhs,A,pivot,b);

  return(status);
}
#endif

/*----------------------------------------------------------------------------*/

int Bsolve(double *A,int hbw, int neq,int *pivot, double *b)

/*----------------------------------------------------------------------------*/
{
  int status,job,nrhs=1,n=3*hbw+1;
  char cjob='N';

  printf("factorisation\n");


#if LAPACKC == 1
  dgbtrf ( neq, neq, hbw, hbw,A, n,pivot,&status);
#elif LAPACKF_ == 1
  dgbtrf_ ( &neq, &neq, &hbw, &hbw,A, &n, pivot, &status);
#elif ATLAS == 1
  printf(" dgbtrf (factorisation) can not be realize using Atlas --> Return Bsolve not realize\n");
#endif

  if(status!=0) return(status);
  printf("solving\n");

#if LAPACKC == 1
  dgbtrs ('N', neq, hbw, hbw, nrhs,A, n,pivot,b, neq,&status);
#elif LAPACKF_ == 1
  dgbtrs_ (&cjob, &neq, &hbw, &hbw, &nrhs,A, &n,pivot,b, &neq,&status);
#elif ATLAS == 1
  printf(" dgbtrs (solve) can not be realize using Atlas --> Return Bsolve not realize\n");
#endif


  return(status);
}

/*----------------------------------------------------------------------------*/

int Psolve(double *A,int hbw, int neq,int *pivot, double *b)

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------

----------------------------------------------------------------------*/
{
  int status,job,nrhs=1,n=hbw+1;
  char cjob='N',cflag='U';

  printf("factorisation\n");

#if LAPACKC == 1
  dpbsv(cflag, neq, hbw, nrhs, A, n, b, neq, &status);
#elif LAPACKF_ == 1
  dpbsv_(&cflag, &neq, &hbw, &nrhs, A,&n, b, &neq, &status);
#elif ATLAS == 1
  printf(" dpbsv  can not be realize using Atlas --> Return Psolve not realize\n");
#endif
 
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int readmatrix(char *matfile,sparse_t *sparse)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,m,n,nx,ny;
  int nitems,status,ncmax=0,nrmax=0,*nc,*nr;
  int *col,*row;
  FILE *in;
  double dum,*matrix;

/*-----------------------------------------------------------------------------
give_psizet_surd.m  mat_cons.dat  matlab_psizet.dat  rhs_cons.dat  */


  in=fopen(matfile,"r");

/*-----------------------------------------------------------------------------
  matrix is given per row, i.e. transpose given by column */
  nx=0;
  ny=0;
  while (!feof(in))
    {
    nitems=fscanf(in,"%d %d %lf",&i,&j,&dum);
    if(nitems!=3) break;
    nx=MAX(i,nx);
    ny=MAX(j,ny);
    }

  rewind(in);
  printf("dimension: nx=%d ny=%d\n",nx,ny);

/*-----------------------------------------------------------------------------
  matrix is given per row, i.e. transpose given by column */

  nc=(int *)malloc(nx*sizeof(int));
  nr=(int *)malloc(ny*sizeof(int));

  for(i=0;i<nx;i++) nc[i]=0;
  for(j=0;j<ny;j++) nr[j]=0;

  while (!feof(in))
    {
    nitems=fscanf(in,"%d %d %lf",&i,&j,&dum);
    if(nitems!=3) break;
    i--;
    j--;
    nc[i]++;
    nr[j]++;
    ncmax=MAX(ncmax,nc[i]);
    nrmax=MAX(nrmax,nr[j]);
    }
  printf("max non-zero row per column=%d\n",nrmax);
  printf("max non-zero column per row=%d\n",ncmax);

  rewind(in);

/*-----------------------------------------------------------------------------
  read sparse matrix A per column*/

  matrix=(double *)malloc(nrmax*ny*sizeof(double));
  row=(int *)malloc(nrmax*ny*sizeof(int));
  col=(int *)malloc(ncmax*nx*sizeof(int));

  for(i=0;i<nx;i++) nc[i]=0;
  for(j=0;j<ny;j++) nr[j]=0;

  printf("read matrix A in sparse format\n");
  while (!feof(in))
    {
    nitems=fscanf(in,"%d %d %lf",&i,&j,&dum);
/*     printf("%d %d %lf\n",i,j,dum); */
    if(nitems!=3) break;
    i--;
    j--;
    if(nr[j]==nrmax) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
    matrix[nrmax*j+nr[j]]=dum;
    row   [nrmax*j+nr[j]]=i;
    col   [ncmax*i+nc[i]]=j;
    nr[j]++;
    nc[i]++;
    }

  fclose(in);

  sparse->nx=nx;
  sparse->ny=ny;
  sparse->ncmax=ncmax;
  sparse->nrmax=nrmax;
  sparse->nc=nc;
  sparse->nr=nr;
  sparse->row=row;
  sparse->col=col;
  sparse->matrix=matrix;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double *readrhs(char *rhsfile,sparse_t sparse)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,m,n,nx,ny;
  int nitems,status,ncmax=0,nrmax=0,*nc,*nr;
  int *col,*row;
  FILE *in;
  double dum,*matrix,*rhs,*vector;

  rhs=(double *)malloc(sparse.ny*sizeof(double));
  vector=(double *)malloc(sparse.nx*sizeof(double));

  for(m=0;m<sparse.ny;m++) rhs[m]=0;

  printf("read rhs file\n");
  in=fopen(rhsfile,"r");

  i=0;
  while (!feof(in))
    {
    nitems=fscanf(in,"%lf",&vector[i]);
    if(nitems!=1) break;
    i++;
    }

  fclose(in);

  printf("compute tA y\n");
  for(j=0;j<sparse.ny;j++) {
    for(m=0;m<sparse.nr[j];m++) {
      i=sparse.row[sparse.nrmax*j+m];
      rhs[j]+=sparse.matrix[sparse.nrmax*j+m]*vector[i];
      }
    }

  free(vector);

  return(rhs);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int product(sparse_t sparse1, sparse_t *sparse2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,m,n,nx,ny;
  int nitems,status,nc1max=0,nr1max=0,nc2max=0,nr2max=0,*nc1,*nr1,*nc2,*nr2;
  int *col1,*row1, *col2,*row2, *colmax,r1,r2;
  FILE *in;
  double dum,*matrix1,*matrix2,*vector,tmp;

/*-----------------------------------------------------------------------------
 */
  nx=sparse1.nx;
  ny=sparse1.ny;
  nc1max=sparse1.ncmax;
  nr1max=sparse1.nrmax;
  nc1=sparse1.nc;
  nr1=sparse1.nr;
  row1=sparse1.row;
  col1=sparse1.col;
  matrix1=sparse1.matrix;


/*-----------------------------------------------------------------------------
  find out the column rightward interaction band*/

  colmax  =(int *)malloc(ny*sizeof(int));
  for(j=0;j<ny;j++) colmax[j]=0;

  for(j=0;j<ny;j++) {
    for(m=0;m<nr1[j];m++) {
      i=row1[nr1max*j+m];
      for(n=0;n<nc1[i];n++) {
        k=col1[nc1max*i+n];
        colmax[j]=MAX(colmax[j],k-j);
        }
      }
    }

/*-----------------------------------------------------------------------------
  compute tAxA  */

  vector=(double *)malloc(nx*sizeof(double));

  nr2=(int *)malloc(ny*sizeof(int));
  nc2=(int *)malloc(ny*sizeof(int));

  for(i=0;i<ny;i++) nc2[i]=0;
  for(j=0;j<ny;j++) nr2[j]=0;

  matrix2=(double *)malloc(2*nr1max*ny*sizeof(double));
  row2   =(int *)malloc(2*nr1max*ny*sizeof(int));
  col2   =(int *)malloc(2*nr1max*ny*sizeof(int));

  for(i=0;i<nx;i++) vector[i]=0;

  printf("compute tA x A\n");
  for(j=0;j<ny;j++) {
    tmp=0;
    for(m=0;m<nr1[j];m++) {
      vector[row1[nr1max*j+m]]=matrix1[nr1max*j+m];
      tmp+=matrix1[nr1max*j+m]*matrix1[nr1max*j+m];
      }
    if(tmp!=0) {
      matrix2[2*nr1max*j+nr2[j]]=tmp; /* (j,j) of tA x A */
      row2[2*nr1max*j+nr2[j]]=j;
      col2[2*nr1max*j+nc2[j]]=j;
      nr2[j]++;
      nc2[j]++;
      }
    if(j%1000==0) printf("%6d over %d\n",j,ny);
    r1=row1[nr1max*j];
    r2=row1[nr1max*j+nr1[j]-1];
/*-----------------------------------------------------------------------------
    final matrix is symmetric !!! */
    for(i=j+1;i<j+colmax[j]+1;i++) {
      if(row1[nr1max*i]>r2) continue;
      if(row1[nr1max*i+nr1[i]-1]<r1) continue;
      tmp=0;
      for(m=0;m<nr1[i];m++) {
        k=row1[nr1max*i+m];
        if(k>r2) break;
        tmp+=matrix1[nr1max*i+m]*vector[k];
        }
      if(tmp!=0) {
        matrix2[2*nr1max*j+nr2[j]]=tmp; /* (i,j) of tA x A */
        row2[2*nr1max*j+nr2[j]]=i;
        col2[2*nr1max*i+nc2[i]]=j;
        nr2[j]++;
        nc2[i]++;
/*-----------------------------------------------------------------------------
        final matrix is symmetric !!! */
        matrix2[2*nr1max*i+nr2[i]]=tmp; /* (j,i) of tA x A */
        row2[2*nr1max*i+nr2[i]]=j;
        col2[2*nr1max*j+nc2[j]]=i;
        nr2[i]++;
        nc2[j]++;
        }
      }
    for(m=0;m<nr1[j];m++) vector[row1[nr1max*j+m]]=0.;
    nr2max=MAX(nr2max,nr2[j]);
    nc2max=MAX(nc2max,nc2[j]);
    }

  sparse2->nx=ny;
  sparse2->ny=ny;
  sparse2->ncmax=nc2max;
  sparse2->nrmax=nr2max;
  sparse2->nc=nc2;
  sparse2->nr=nr2;
  sparse2->row=row2;
  sparse2->col=col2;
  sparse2->matrix=matrix2;

  free(vector);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int solve01(char *matfile, char *rhsfile, char *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,m,n,nx,ny,bw,hbw,previous;
  int nitems,status,ncmax=0,nrmax=0,*nc,*nr;
  int ml,mu;
  int *row_skyline[2],*col_skyline[2];
  int *col,*row;
  int *colmax;
  int r1,r2;
  FILE *in;
  double *mat,*rhs,dum,*matrix,*vector;
  sparse_t sparse,tmp;

/*-----------------------------------------------------------------------------
  read matrix and initialise row and column index arrays*/

  readmatrix(matfile,&tmp);
 
  printf("max non-zero row per column=%d\n",tmp.nrmax);
  printf("max non-zero column per row=%d\n",tmp.ncmax);

  rhs=readrhs(rhsfile,tmp);

  product(tmp,&sparse);
  
  nx=sparse.nx;
  ny=sparse.ny;
  ncmax=sparse.ncmax;
  nrmax=sparse.nrmax;
  nc=sparse.nc;
  nr=sparse.nr;
  row=sparse.row;
  col=sparse.col;
  matrix=sparse.matrix;

/*-----------------------------------------------------------------------------
  lower band width */
  hbw=0;
  for(j=0;j<ny;j++) {
    for(m=0;m<nr[j];m++) {
      hbw=(int)( MAX(hbw,fabs( 1.0* j-row[2*tmp.nrmax*j+m]))  ); //le 1.0 est necessaire pour passer l'expression en double
      }
    }

  bw=2*hbw+1;
  printf("bandwidth: %d %d %d\n",hbw,bw,nrmax);

/*-----------------------------------------------------------------------------
  lower band width */
  ml=0;
  for(j=0;j<ny;j++) {
    for(m=0;m<nr[j];m++) {
      ml=MAX(ml,row[2*tmp.nrmax*j+m]-j);
      }
    }

/*-----------------------------------------------------------------------------
  upper band width */
  mu=0;
  for(i=0;i<ny;i++) {
    for(m=0;m<nc[i];m++) {
      mu=MAX(mu,col[2*tmp.nrmax*i+m]-i);
      }
    }
  printf("bandwidth: ml=%d mu=%d\n",ml,mu);

  mat=(double *)malloc((3*hbw+1)*ny*sizeof(double));

  for(m=0;m<(3*hbw+1)*ny;m++) mat[m]=0;

  for(j=0;j<ny;j++) {
    for(m=0;m<nr[j];m++) {
      i=row[2*tmp.nrmax*j+m];
      mat[(3*hbw+1)*j+i-j+bw-1]=matrix[2*tmp.nrmax*j+m];
      }
    }

/*-----------------------------------------------------------------------------
  linear system solving */
  printf("solve system\n");
  status=Bsolve(mat, hbw, ny,nr,rhs);
  printf("status=%d\n",status);

  printf("save solution file\n");
  in=fopen(output,"w");
  for(j=0;j<ny;j++) {
    fprintf(in,"%lf\n",rhs[j]);
    }

  fclose(in);
  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int solve02(char *matfile, char *rhsfile, char *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,m,n,nx,ny,bw,hbw,previous;
  int nitems,status,ncmax=0,nrmax=0,*nc,*nr;
  int ml,mu;
  int *pivot;
  int *col,*row;
  int *colmax;
  int r1,r2;
  FILE *in;
  double *mat,*rhs,dum,*matrix,*vector;
  sparse_t sparse,tmp;

/*-----------------------------------------------------------------------------
  read matrix and initialise row and column index arrays*/

  readmatrix(matfile,&tmp);
 
  printf("max non-zero row per column=%d\n",tmp.nrmax);
  rhs=readrhs(rhsfile,tmp);

  printf("max non-zero column per row=%d\n",tmp.ncmax);

  product(tmp,&sparse);
  
  nx=sparse.nx;
  ny=sparse.ny;
  ncmax=sparse.ncmax;
  nrmax=sparse.nrmax;
  nc=sparse.nc;
  nr=sparse.nr;
  row=sparse.row;
  col=sparse.col;
  matrix=sparse.matrix;

/*-----------------------------------------------------------------------------
  lower band width */
  hbw=0;
  for(j=0;j<ny;j++) {
    for(m=0;m<nr[j];m++) {
      hbw=(int)( MAX(hbw,fabs(1.00 * j-row[2*tmp.nrmax*j+m])) );
      }
    }

  bw=2*hbw+1;
  printf("bandwidth: %d %d %d\n",hbw,bw,nrmax);

/*-----------------------------------------------------------------------------
  lower band width */
  ml=0;
  for(j=0;j<ny;j++) {
    for(m=0;m<nr[j];m++) {
      ml=MAX(ml,row[2*tmp.nrmax*j+m]-j);
      }
    }

/*-----------------------------------------------------------------------------
  upper band width */
  mu=0;
  for(i=0;i<ny;i++) {
    for(m=0;m<nc[i];m++) {
      mu=MAX(mu,col[2*tmp.nrmax*i+m]-i);
      }
    }
  printf("bandwidth: ml=%d mu=%d\n",ml,mu);

  mat=(double *)malloc((hbw+1)*ny*sizeof(double));
  pivot=(int *)malloc(ny*sizeof(int));

  for(m=0;m<(hbw+1)*ny;m++) mat[m]=0;

  for(j=0;j<ny;j++) {
    for(m=0;m<nr[j];m++) {
      i=row[2*tmp.nrmax*j+m];
      if(i<=j)
        mat[(hbw+1)*j+i-j+hbw]=matrix[2*tmp.nrmax*j+m];
      }
    }

/*-----------------------------------------------------------------------------
  linear system solving */
  printf("solve system\n");
  status=Psolve(mat, hbw, ny, pivot, rhs);
  printf("status=%d\n",status);

  printf("save solution file\n");
  in=fopen(output,"w");
  for(j=0;j<ny;j++) {
    fprintf(in,"%lf\n",rhs[j]);
    }

  fclose(in);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int solve03(char *matfile, char *rhsfile, char *output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,m,n,nx,ny,bw,hbw,previous;
  int nitems,status,ncmax=0,nrmax=0,*nc,*nr,*nc3,*nr3,nc3max=0,nr3max=0;
  int ml,mu;
  int *row_skyline[2],*col_skyline[2];
  int *col2,*row2,*col3,*row3;
  int *colmax;
  int r1,r2;
  FILE *in;
  double *mat,*rhs,dum,*sparse1,*sparse2,*sparse3,*vector,tmp;

/*-----------------------------------------------------------------------------
give_psizet_surd.m  mat_cons.dat  matlab_psizet.dat  rhs_cons.dat  */


  in=fopen("../symphonie/mat_cons.dat","r");

/*-----------------------------------------------------------------------------
  matrix is given per row, i.e. transpose given by column */
  nx=0;
  ny=0;
  while (!feof(in))
    {
    nitems=fscanf(in,"%d %d %lf",&i,&j,&dum);
    if(nitems!=3) break;
    nx=MAX(i,nx);
    ny=MAX(j,ny);
    }

  rewind(in);
  printf("dimension: nx=%d ny=%d\n",nx,ny);

/*-----------------------------------------------------------------------------
  matrix is given per row, i.e. transpose given by column */

  nc=(int *)malloc(nx*sizeof(int));
  nr=(int *)malloc(ny*sizeof(int));

  nc3=(int *)malloc(ny*sizeof(int));
  nr3=(int *)malloc(ny*sizeof(int));

  for(i=0;i<nx;i++) nc[i]=0;
  for(j=0;j<ny;j++) nr[j]=0;
  while (!feof(in))
    {
    nitems=fscanf(in,"%d %d %lf",&i,&j,&dum);
    if(nitems!=3) break;
    i--;
    j--;
    nc[i]++;
    nr[j]++;
    ncmax=MAX(ncmax,nc[i]);
    nrmax=MAX(nrmax,nr[j]);
    }
  printf("max non-zero row per column=%d\n",nrmax);
  printf("max non-zero column per row=%d\n",ncmax);

  rewind(in);

/*-----------------------------------------------------------------------------
  read sparse matrix A per column*/

  sparse2=(double *)malloc(nrmax*ny*sizeof(double));
  row2=(int *)malloc(nrmax*ny*sizeof(int));
  col2=(int *)malloc(ncmax*nx*sizeof(int));

  for(i=0;i<nx;i++) nc[i]=0;
  for(j=0;j<ny;j++) nr[j]=0;

  printf("read matrix A in sparse format\n");
  while (!feof(in))
    {
    nitems=fscanf(in,"%d %d %lf",&i,&j,&dum);
/*     printf("%d %d %lf\n",i,j,dum); */
    if(nitems!=3) break;
    i--;
    j--;
    if(nr[j]==nrmax) {__ERR_BASE_LINE__("exiting\n");exit(-1);}
    sparse2[nrmax*j+nr[j]]=dum;
    row2   [nrmax*j+nr[j]]=i;
    col2   [ncmax*i+nc[i]]=j;
    nr[j]++;
    nc[i]++;
    }

  fclose(in);

/*-----------------------------------------------------------------------------
  find out the column rightward interaction band*/

  colmax  =(int *)malloc(ny*sizeof(int));
  for(j=0;j<ny;j++) colmax[j]=0;

  for(j=0;j<ny;j++) {
    for(m=0;m<nr[j];m++) {
      i=row2[nrmax*j+m];
      for(n=0;n<nc[i];n++) {
        k=col2[ncmax*i+n];
        colmax[j]=MAX(colmax[j],k-j);
        }
      }
    }

/*-----------------------------------------------------------------------------
  compute tAxA  */

  vector=(double *)malloc(nx*sizeof(double));
  sparse3=(double *)malloc(2*nrmax*ny*sizeof(double));
  row3   =(int *)malloc(2*nrmax*ny*sizeof(int));
  col3   =(int *)malloc(2*nrmax*ny*sizeof(int));

  for(j=0;j<ny;j++) nr3[j]=0;
  for(i=0;i<ny;i++) nc3[i]=0;
  for(i=0;i<nx;i++) vector[i]=0;

  printf("compute tA x A\n");
  for(j=0;j<ny;j++) {
    tmp=0;
    for(m=0;m<nr[j];m++) {
      vector[row2[nrmax*j+m]]=sparse2[nrmax*j+m];
      tmp+=sparse2[nrmax*j+m]*sparse2[nrmax*j+m];
      }
    if(tmp!=0) {
      sparse3[2*nrmax*j+nr3[j]]=tmp; /* (j,j) of tA x A */
      row3[2*nrmax*j+nr3[j]]=j;
      col3[2*nrmax*j+nc3[j]]=j;
      nr3[j]++;
      nc3[j]++;
      }
    if(j%1000==0) printf("%6d over %d\n",j,ny);
    r1=row2[nrmax*j];
    r2=row2[nrmax*j+nr[j]-1];
/*-----------------------------------------------------------------------------
    final matrix is symmetric !!! */
    for(i=j+1;i<j+colmax[j]+1;i++) {
      if(row2[nrmax*i]>r2) continue;
      if(row2[nrmax*i+nr[i]-1]<r1) continue;
      tmp=0;
      for(m=0;m<nr[i];m++) {
        k=row2[nrmax*i+m];
        if(k>r2) break;
        tmp+=sparse2[nrmax*i+m]*vector[k];
        }
      if(tmp!=0) {
        sparse3[2*nrmax*j+nr3[j]]=tmp; /* (i,j) of tA x A */
        row3[2*nrmax*j+nr3[j]]=i;
        col3[2*nrmax*i+nc3[i]]=j;
        nr3[j]++;
        nc3[i]++;
/*-----------------------------------------------------------------------------
        final matrix is symmetric !!! */
        sparse3[2*nrmax*i+nr3[i]]=tmp; /* (j,i) of tA x A */
        row3[2*nrmax*i+nr3[i]]=j;
        col3[2*nrmax*j+nc3[j]]=i;
        nr3[i]++;
        nc3[j]++;
        }
      }
    for(m=0;m<nr[j];m++) vector[row2[nrmax*j+m]]=0.;
    nr3max=MAX(nr3max,nr3[j]);
    nc3max=MAX(nc3max,nc3[j]);
    }

  printf("max non-zero row per column=%d\n",nr3max);
  printf("max non-zero column per row=%d\n",nc3max);

/*
  row_skyline[0]=malloc(nx*sizeof(int));
  row_skyline[1]=malloc(nx*sizeof(int));

  col_skyline[0]=malloc(ny*sizeof(int));
  col_skyline[1]=malloc(ny*sizeof(int));
*/
/*-----------------------------------------------------------------------------
  lower band width */
  hbw=0;
  for(j=0;j<ny;j++) {
    for(m=0;m<nr3[j];m++) {
      hbw=(int)( MAX(hbw,fabs(1.0 * j-row3[2*nrmax*j+m])) );
      }
    }

  bw=2*hbw+1;
  printf("bandwidth: %d %d %d\n",hbw,bw,nr3max);

/*-----------------------------------------------------------------------------
  lower band width */
  ml=0;
  for(j=0;j<ny;j++) {
    for(m=0;m<nr3[j];m++) {
      ml=MAX(ml,row3[2*nrmax*j+m]-j);
      }
    }

/*-----------------------------------------------------------------------------
  upper band width */
  mu=0;
  for(i=0;i<ny;i++) {
    for(m=0;m<nc3[i];m++) {
      mu=MAX(mu,col3[2*nrmax*i+m]-i);
      }
    }
  printf("bandwidth: ml=%d mu=%d\n",ml,mu);

  mat=(double *)malloc((3*hbw+1)*ny*sizeof(double));
  rhs=(double *)malloc(ny*sizeof(double));

  for(m=0;m<(3*hbw+1)*ny;m++) mat[m]=0;
  for(m=0;m<ny;m++) rhs[m]=0;

  for(j=0;j<ny;j++) {
    for(m=0;m<nr3[j];m++) {
      i=row3[2*nrmax*j+m];
      mat[(3*hbw+1)*j+i-j+bw-1]=sparse3[2*nrmax*j+m];
      }
    }

  printf("read rhs file\n");
  in=fopen("../symphonie/rhs_cons.dat","r");

  i=0;
  while (!feof(in))
    {
    nitems=fscanf(in,"%lf",&vector[i]);
    if(nitems!=1) break;
    i++;
    }

  fclose(in);

  printf("compute tA y\n");
  for(j=0;j<ny;j++) {
    for(m=0;m<nr[j];m++) {
      i=row2[nrmax*j+m];
      rhs[j]+=sparse2[nrmax*j+m]*vector[i];
      }
    }

/*-----------------------------------------------------------------------------
  linear system solving */
  printf("solve system\n");
  status=Bsolve(mat, hbw, ny,nr,rhs);
  printf("status=%d\n",status);

  printf("save solution file\n");
  in=fopen(output,"w");
  for(j=0;j<ny;j++) {
    fprintf(in,"%lf\n",rhs[j]);
    }

  fclose(in);
  printf("save solution file\n");
  in=fopen("../symphonie/tutu.dat","w");
  for(j=0;j<ny;j++) {
    for(m=0;m<nr3[j];m++) {
      i=row3[2*nrmax*j+m];
      fprintf(in,"%d %d %lf\n",i,j,sparse3[2*nrmax*j+m]);
      
      }
    }

  fclose(in);
  return(0);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status;

  size_t size;
  FILE *in,*out;
  char *keyword,*zone,*s;

  char *rootname=NULL,*output=NULL,*input=NULL,*pathname=NULL,*bathymetry=NULL;
  char *bathy=NULL,*notebook=NULL,*rhs=NULL,*matrice=NULL;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
        case 'm' :
          matrice= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rhs= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        __OUT_BASE_LINE__("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
      free(keyword);
    }

  input=(char* )malloc(1024);
  output=(char* )malloc(1024);

  if(matrice ==NULL) matrice=strdup("../symphonie/mat_cons.dat");
  if(rhs ==NULL)     rhs=strdup("../symphonie/rhs_cons.dat");
  if(output ==NULL)  output=strdup("../symphonie/toto.dat");
/*
  solve01(matrice,rhs,output);
*/
  solve02(matrice,rhs,output);
/*
  solve03(matrice,rhs,output);
*/
  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
  
}
