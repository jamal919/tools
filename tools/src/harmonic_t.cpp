
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief harmonic_t methods that use solvers
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string>

#include "tides.h"

#include "admittance-mts.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_t::factorize(int destructive)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nw,status;
  
  nw=spectrum.n;
  if(nw<=0 || A==0) return -1;
  
  neq=2*nw;
  pivot = new int[neq];
  
  if(destructive==1) {
    M = A;
    }
  else {
    M = new double[neq * neq];
    for(int k = 0; k < neq * neq; k++) {
      M[k] = A[k];
      }
    }
  status=poc_getrf(neq, M, pivot);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_t::solve(double *b, int nrhs, int destructive)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int job, status=0;
  
  for(int step = 0; step < 1; step++) {
    job = 0;
    int nrhs=1;
    status=poc_getrs(neq, nrhs, M, pivot, b);
    if(status!=0)
      return status;
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_t::error(double *residuals, int n, double* & error)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  nrhs=1;
  
  double s2 = 0;
  for(size_t i = 0; i < n; i++){
    s2 += residuals[i] * residuals[i];
    }
/* *----------------------------------------------------------------------------
  unbiased estimator of the residual variance */
  const int
    unbiasing=n - neq - 1;
  
  const double
    rus2=sqrt(s2/unbiasing);
  
  if(not isfinite(rus2))
    TRAP_ERR_RETURN(-1,1,"ERROR : unbiased estimator of the residual variance is %g=sqrt(%g/%d=(%d-%d-1))\n",
      rus2,s2,unbiasing,n,neq);

  deletep(&error);
  error = new double[neq];
  
  double *mu    = new double[neq];
  double *bck   = new double[neq];
  for(size_t i = 0; i < neq; i++){
    double c = 0.;
    for(size_t j = 0; j < neq; j++) mu[j]=0.0;
/* *----------------------------------------------------------------------------
    commentaires ? */
    mu[i] = 1.;
    memcpy(bck, mu, neq * sizeof(double));
/* *----------------------------------------------------------------------------
    modified normal equations */
    status=poc_getrs(neq, nrhs, M, pivot, mu);
    if(status!=0)return status;
    for(size_t j = 0; j < neq; j++){
      c += bck[j] * mu[j];
      }
/* *----------------------------------------------------------------------------
    error on i-th element estimates */    
    if(!isfinite(c)){
      printf("alert %d\n",i);
      }

    error[i]=c*rus2;
    }
  delete[] mu;
  delete[] bck;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void matrix_member(const double *A, const spectrum_t & s, int k, int l, double *rere, double *imim, double *reim, double *imre)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const int neq = 2 * s.n;
  
  const int i=2 * l * neq  + 2 * k;
  
  *rere=A[i];
  
  *imim=A[i+neq+1];
  
  if(reim!=0)
    *reim=A[i+neq];
  
  if(imre!=0)
    *imre=A[i+1];
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void matrix_auto_correlation(const double *A, const spectrum_t & s, int k, double *re, double *im)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double rere,imim;
  
  matrix_member(A,s,k,k,&rere,&imim);
  
  *re=sqrt(rere);
  *im=sqrt(imim);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double matrix_correlation(const double *A, const spectrum_t & s, int k, int l, double rk, double ik)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double rl,il;
  matrix_auto_correlation(A,s,l,&rl,&il);
  
  double rkrl,ikil,rkil,ikrl;
  matrix_member(A,s,k,l,&rkrl,&ikil,&rkil,&ikrl);
  
  double m=0.;
  updatemax(&m,abs(rkrl/(rk*rl)));
  updatemax(&m,abs(ikil/(ik*il)));
  updatemax(&m,abs(rkil/(rk*il)));
  updatemax(&m,abs(ikrl/(ik*rl)));
  
  if(m>1.1){
    STDERR_BASE_LINE("k %10.4g:%10.4g %10.4g:%10.4g\n",rk,rk*rk,ik,ik*ik);
    STDERR_BASE_LINE("l %10.4g:%10.4g %10.4g:%10.4g\n",rl,rl*rl,il,il*il);
    STDERR_BASE_LINE("%10.4g/(%10.4g*%10.4g)=%10.4g\n",rkrl,rk,rl,rkrl/(rk*rl));
    STDERR_BASE_LINE("%10.4g/(%10.4g*%10.4g)=%10.4g\n",ikil,ik,il,ikil/(ik*il));
    STDERR_BASE_LINE("%10.4g/(%10.4g*%10.4g)=%10.4g\n",rkil,rk,il,rkil/(rk*il));
    STDERR_BASE_LINE("%10.4g/(%10.4g*%10.4g)=%10.4g\n",ikrl,ik,rl,ikrl/(ik*rl));
    TRAP_ERR_EXIT(ENOEXEC,"programming error : %10.4g\n",m);
    }
  
  return m;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T>  void save_vector_template(FILE *f, const char *title, const char *newLine, const char *newCell, const char *format, const T *values, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(values==0) return;
  bool skip=n<0;
  if(skip)
    n=-n;
  
  fprintf(f,"%s%s",newLine,title);
  
  for(int k=0;k<n;k++){
    fprintf(f,format,newCell,values[k]);
    if(skip)
      fprintf(f,"%s",newCell);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void save_vector(FILE *f, const char *title, const char *newLine, const char *newCell, const int *values, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  save_vector_template(f,title,newLine,newCell,"%s%d",values,n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void save_vector(FILE *f, const char *title, const char *newLine, const char *newCell, const double *values, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  save_vector_template(f,title,newLine,newCell,"%s%.3g",values,n);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void save_full_matrix(const string & path, int *index, int tag, const spectrum_t & s, const double *A, const int *keep, const int *deduce, double *rhs, double *zr, double *zi, double *error)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// save matrix, RHS, solution and many optional more to CSV or HTML
/**
When the path ends with .htm(l), select HTML format. CSV format otherwise.
\note this can show correlations >1 after admittance!
*/
/*----------------------------------------------------------------------------*/
{
  const int neq = 2 * s.n;
  size_t k,l,i;
  
  FILE *f;
  string mode="w";
#define TABLE_FORMAT_UNDEFINED  0
#define TABLE_FORMAT_HTML 1
#define TABLE_FORMAT_CSV  2
  int format=TABLE_FORMAT_UNDEFINED;
  
  if(strrncasecmp(path,".html")==0 || strrncasecmp(path,".htm")==0)
    format=TABLE_FORMAT_HTML;
  else // if(strrncasecmp(path,".csv")==0)
    format=TABLE_FORMAT_CSV;
  
  /* TODO: path mutex instead of pragma critical */
  #pragma omp critical(save_full_matrix_mutex)
  {
  if(index!=0){
    if(*index>0)
      mode[0]='a';
    
    (*index)++;
    }
  }
  
  const int max_index=1000;
  if(*index>max_index) TRAP_ERR_RETURN(,1,"More than %d matrices already saved to "+path+". Skipping %d th.\n",max_index,*index);
  
  #pragma omp critical(save_full_matrix_mutex)
  {
  switch(mode[0]){
  case 'a':
    printf("Appending matrix to "+path+"\n");
    break;
  case 'w':
    printf("Overwriting "+path+" with matrix.\n");
    }
  f=fopen(path.c_str(),mode.c_str());
  
/*------------------------------------------------------------------------------
  header */
  const char *newLine,*newCell;
  switch(format){
  case TABLE_FORMAT_HTML:
    fprintf(f,"\n<p>Station No%5d",tag);
    fprintf(f,"<table>");
    
    newLine="\n<tr><td>";
    newCell="<td>";
    break;
  case TABLE_FORMAT_CSV:
    if(mode[0]=='w')
      fprintf(f,"Use format code:,[MAGENTA][<-0.3]#0.##;[RED][>0.3]#0.##;[BLACK]#0.##");
    fprintf(f,"\nStation No%5d",tag);
    
    newLine="\n";
    newCell=",";
    break;
    }
  
  fprintf(f,"%swave",newLine);
  for(k = 0; k < neq; k++) {
    fprintf(f,"%s%s",newCell,s.waves[k/2].name);
    }
  
/*------------------------------------------------------------------------------
  matrix */
  double *diagonal=new double[neq];
  for(k = 0; k < neq; k++) {
    diagonal[k]=A[k*neq+k];
    }
  
  for(k = 0,i=0; k < neq; k++) {
    fprintf(f,"%s%s",newLine,s.waves[k/2].name);
    for(l = 0; l < neq; l++,i++) {
      fprintf(f,"%s%.3g",newCell,A[i]/sqrt(diagonal[k]*diagonal[l]));
      }
    }
  
/*------------------------------------------------------------------------------
  vectors */
  save_vector(f,"diagonal",newLine,newCell,diagonal,neq);
  deletep(&diagonal);
  
  save_vector(f,"keep",newLine,newCell,keep,-s.n);
  save_vector(f,"deduce",newLine,newCell,deduce,-s.n);
  
  save_vector(f,"RHS",newLine,newCell,rhs,neq);
  
  if(zr!=0 && zi!=0){
    double solution[neq];
    for(k=0,i=0; k < neq; k+=2,i++){
      solution[k]=zr[i];
      solution[k+1]=zi[i];
      }
    save_vector(f,"solution",newLine,newCell,solution,neq);
    }
  
  save_vector(f,"error",newLine,newCell,error,neq);
  
/*------------------------------------------------------------------------------
  tail */
  switch(format){
  case TABLE_FORMAT_HTML:
    fprintf(f,"\n</table>\n");
    break;
    }
  
  fclose(f);
  }
}
