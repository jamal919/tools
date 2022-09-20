
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief definition of VERY OLD utility functions
Other utility functions are defined in functions.cpp
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <iostream>
#include <stdlib.h>
#include <string.h>

#include "poc-assertions.h"
#include "tools-define.h"

#include "rutin.h"

/*----------------------------------------------------------------------------*/

void gmerror(const char *error_text)

/*----------------------------------------------------------------------------*/
{
/*  void __ERR_BASE_LINE__("exiting\n");exit(int);     */

  fprintf(stderr,"Execution terminated by checks in program...\n");
  fprintf(stderr,"%s\n",error_text);
  __ERR_BASE_LINE__("...now exiting to system...\n");
  exit(1);
}

/*----------------------------------------------------------------------------*/

void gmerror(string error_text)

/*----------------------------------------------------------------------------*/
{
  cerr << "Execution terminated by checks in program..." << endl;
  cerr << error_text << endl;
  cerr << "...now exiting to system..." << endl;
  __ERR_BASE_LINE__("exiting\n");exit(1);
}


/*----------------------------------------------------------------------------*/

void bread_dg(char *addr, size_t size, size_t n, FILE *file)

/*----------------------------------------------------------------------------*/
{
  if(fread(addr,size,n,file)!= n){
        gmerror("read error bread_dg");
  }
}
 
/*----------------------------------------------------------------------------*/

void bwrite_dg(char *addr, size_t size, size_t n, FILE *file)

/*----------------------------------------------------------------------------*/
{
  if(fwrite(addr,size,n,file)!= n){
        gmerror("write error bwrite_dg");
  }
}
 
  
/*----------------------------------------------------------------------------*/

  void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch)

/*----------------------------------------------------------------------------*/
{
  int i;
 
  for(i=nrh;i>=nrl;i--) free( (m[i]+ncl));
  free( (m+nrl));
}
 
/*----------------------------------------------------------------------------*/

  void free_ivector(int *v,int nl,int nh)

/*----------------------------------------------------------------------------*/
{
  free( (v+nl));
}
 
/*----------------------------------------------------------------------------*/

  void free_smatrix(float **m,int nrl,int nrh,int ncl,int nch)

/*----------------------------------------------------------------------------*/
{
  int i;
 
  for(i=nrh;i>=nrl;i--) free( (m[i]+ncl));
  free( (m+nrl));
}
 
/*----------------------------------------------------------------------------*/

  void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch)

/*----------------------------------------------------------------------------*/
{
  int i;
 
  for(i=nrh;i>=nrl;i--) free( (m[i]+ncl));
  free( (m+nrl));
}
/*----------------------------------------------------------------------------*/

  void free_svector(float *v,int nl,int nh)

/*----------------------------------------------------------------------------*/
{
  free( (v+nl));
}

/*----------------------------------------------------------------------------*/

  void free_dvector(double *v,int nl,int nh)

/*----------------------------------------------------------------------------*/
{
  free( (v+nl));
}

 
/*----------------------------------------------------------------------------*/

void gmplot(float yy[],unsigned short n,double x1,double dx,FILE *file)

/*----------------------------------------------------------------------------*/
{
  int i,j,x_axis;
  int RANGE = 64;   /* adjust RANGE to get maximum plot on page     */
  float scale;		/* RANGE = 64 for 80 column printer and screen  */
  float x = x1;
  float ymax = -1.e35;
  float ymin = 1.e35;
  float *y=NULL;
 
  y = yy - 1;
  for(i=1; i<=n; i++){
        ymax = MAX(ymax,y[i]);
        ymin = MIN(ymin,y[i]);
  }
  if(ymax == ymin){
        fprintf(file,"\n\n\n********* array all == %f ******\n",ymax);
        return;
  }
  fprintf(file,"\f___X__]");
  for(i=1; i<=RANGE; i++)fprintf(file,"_");
  fprintf(file,"[___Y___\n");
  scale = (RANGE-1)/(ymax-ymin);
  x_axis = RANGE +1;
  if(ymax*ymin < 0.0) x_axis = (int)( fabs(ymin)*scale+1.4999 );
  for(i=1; i<=n; i++){
        fprintf(file," %6.2f]",x);
        for(j=1; j<=RANGE; j++){
                if(j == (int)(scale*(y[i]-ymin)+1.4999)){
                        fprintf(file,"*");
                }else if(j == x_axis){
                        fprintf(file,"|");
                }else{
                        fprintf(file," ");
                }
        }
        fprintf(file,"[%7.2f\n",y[i]);
        x += dx;
  }
}
 
/*----------------------------------------------------------------------------*/

int **imatrix(int nrl,int nrh,int ncl,int nch)

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------

  nrl,nrh: row LOW index,    row HIGH index
  ncl,nch: column LOW index, column HIGH index
  
  allocate a [nrl - nrh, ncl - nch] integer matrix

----------------------------------------------------------------------*/
{
  int i,**matrix=NULL;
 
  exitIfNull(
    matrix=(int **)malloc((size_t) (nrh-nrl+1)*sizeof(int*))
    );
  if (!matrix) gmerror("allocation failure 1 in imatrix()");
  matrix -= nrl;
 
  for(i=nrl;i<=nrh;i++) {
    exitIfNull(
      matrix[i]=(int *)malloc((size_t) (nch-ncl+1)*sizeof(int))
      );
    if (!matrix[i]) gmerror("allocation failure 2 in imatrix()");
    matrix[i] -= ncl;
  }
  return matrix;
}
 
/*----------------------------------------------------------------------------*/

int *ivector(int nl,int nh)

/*----------------------------------------------------------------------------*/
{
  int *v=NULL;
 
  exitIfNull(
    v=(int *)malloc((size_t) (nh-nl+1)*sizeof(int))
    );
  if (!v) gmerror("allocation failure in ivector()");
  return v-nl;
}
 
/*----------------------------------------------------------------------------*/

float **smatrix(int nrl,int nrh,int ncl,int nch)

/*----------------------------------------------------------------------------*/
{
  int i;
  float **m=NULL;
 
  exitIfNull(
    m=(float **) malloc((size_t) (nrh-nrl+1)*sizeof(float*))
    );
  if (!m) gmerror("allocation failure 1 in matrix()");
  m -= nrl;
 
  for(i=nrl;i<=nrh;i++) {
        exitIfNull(
          m[i]=(float *) malloc((size_t) (nch-ncl+1)*sizeof(float))
          );
        if (!m[i]) {
          printf("i, nch-ncl+1: %d %d \n",i, nch-ncl+1);
          gmerror("allocation failure 2 in matrix()");
          }
        m[i] -= ncl;
  }
  return m;
}
 
/*----------------------------------------------------------------------------*/

float *svector(int nl,int nh)

/*----------------------------------------------------------------------------*/
{
  float *v=NULL;
 
  exitIfNull(
    v=(float *)malloc((size_t) (nh-nl+1)*sizeof(float))
    );
  if (!v) gmerror("allocation failure in vector()");
  return v-nl;
}

/*----------------------------------------------------------------------------*/

double **dmatrix(int nrl,int nrh,int ncl,int nch)

/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------

Warning: the returned address m is shifted to the left by a quantity of nrl

----------------------------------------------------------------------*/
{
  int i;
  double **m=NULL;

  exitIfNull(
    m=(double **) malloc((size_t) (nrh-nrl+1)*sizeof(double*))
    );
  if (!m) gmerror("allocation failure 1 in matrix()");
  m -= nrl;
 
  for(i=nrl;i<=nrh;i++) {
        exitIfNull(
          m[i]=(double *) malloc((size_t) (nch-ncl+1)*sizeof(double))
          );
        if (!m[i]) {
          printf("i, nch-ncl+1: %d %d \n",i, nch-ncl+1);
          gmerror("allocation failure 2 in matrix()");
          }
        m[i] -= ncl;
  }
  return m;
}
 
/*----------------------------------------------------------------------------*/

double *dvector(int nl,int nh)

/*----------------------------------------------------------------------------*/
{
  double *v=NULL;

  exitIfNull(
    v=(double *)malloc((size_t) (nh-nl+1)*sizeof(double))
    );
  if (!v) gmerror("allocation failure in vector()");
  return v-nl;
}
 

