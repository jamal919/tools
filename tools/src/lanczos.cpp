
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "tools-structures.h"

#include "tools-define.h"

#include "filter.h"
#include "lanczos.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void lanczos1D_init(double f_c, int l_f, double *weight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc



cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i;
  double q,r,shift,t1,t2;
  double  pi = M_PI;
  
  aset(weight,2*l_f+1,0.);
 
  weight[l_f]=2.*f_c;

  for(i=0;i<l_f;i++) {
    shift=l_f-i;
    t1=shift*f_c;
    r=sin(2.*pi*t1)/(pi*shift);
    t2=shift/l_f;
    q=pi*t2;
    weight[i]=r*sin(q)/q;
    }
  
  for(i=l_f+1;i<2*l_f+1;i++)
    weight[i]=weight[2*l_f-i];
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> void lanczos1D_init_template(T step, int nval, T scale, T **weight, int *nweight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  step  : series interval
  
  nval  : max length of filter
  
  scale : cut-off frequency

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i;
  int n;
  T x,px,y;
  
  n=int(round(scale/step*2));
  n=min(nval,n);
  
  *weight=aset(n,(T) 0.);
  if( *weight == NULL ) n=-1;
  
  for(i=0;i<n;i++) {
    x=i*step/scale/2;
    x=fabs(x);
    if (x <= 1.) {
      if(i==0)
        y=1.;
      else{
        px=x*M_PI;
        y=sin(2*px)*sin(px)/square(px)/2;
        }
      (*weight)[i]=y;
      }
    }
  *nweight=n;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void lanczos1D_init(float step, int nval, float scale, float **weight, int *nweight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  lanczos1D_init_template(step, nval, scale, weight, nweight);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void lanczos1D_init(double step, int nval, double scale, double **weight, int *nweight)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  lanczos1D_init_template(step, nval, scale, weight, nweight);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void lanczos1D_BF(int n, double *h, double mask, double *weight, int m, double *hf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j;
  double sum_w,sum_wh,z;

  for(i=0;i<n;i++) {
    sum_w=0.;
    sum_wh=0.;
    for(j=MAX(i-m,0);j<MIN(i+m,n-1);j++) {
      z=h[j];
      if(z != mask) {
        sum_w=sum_w+weight[m+i-j];
        sum_wh=sum_wh+weight[m+i-j]*z;
        }
      }
    if(sum_w != 0.0)
      hf[i]=sum_wh/sum_w;
    else
      hf[i]=mask;
    }
  for(j=0;j<n;j++) if(h[j] == mask) hf[j]=mask;
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int lanczos1D(int n, double dx, double l_c, double *h, double mask,double *lf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------
   n   : nbre de points dans la serie h                 
   dx  : pas d'echantillonnage de h                     
   l_c : periode de coupure dans la meme unite que dx   
   lf  : signal filtr? (on garde les basses frequences) 
-----------------------------------------------------------------------*/
{
  int  m;
  double *weight=NULL,f_c;
    
  m=MIN((int)(l_c/dx)+1,n);

  f_c=1./((double)(m));
  
  m=4*m;
  
  exitIfNull(
    weight=new double[2*m+1]
    );

  lanczos1D_init(f_c,m,weight);

  lanczos1D_BF(n,h,mask,weight,m,lf);

  delete[] weight;

  return(0);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int lanczos1D(int n, double dx, double *l_c, double *h, double mask,double *lf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------
   n   : nbre de points dans la serie h                 
   dx  : pas d'echantillonnage de h                     
   l_c : periode de coupure dans la meme unite que dx   
   lf  : signal filtr? (on garde les basses frequences) 
-----------------------------------------------------------------------*/
{
  
  for(int i=0;i<n;i++) {
    double *weight=NULL;
    int m=MIN((l_c[i]/dx)+1,n);
    m=4*m;
    weight=new double[2*m+1];
    double sum_w=0.;
    double sum_wh=0.;
    for(int j=MAX(i-m,0);j<MIN(i+m,n-1);j++) {
      double z=h[j];
      if(z != mask) {
        sum_w=sum_w+weight[m+i-j];
        sum_wh=sum_wh+weight[m+i-j]*z;
        }
      }
    if(sum_w != 0.0)
      lf[i]=sum_wh/sum_w;
    else
      lf[i]=mask;
    delete[] weight;
    }
    
  for(int j=0;j<n;j++) if(h[j] == mask) lf[j]=mask;
  
  return(0);

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int Lanczos1D_weights_template(double L, double dx, int a, T * & weight, int & m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  
  http://en.wikipedia.org/wiki/Lanczos_resampling
  
  K = a * sin(pi*x) * sin(pi * x/ a) / (pi *x)?
  
  x is adimensioned distance x=X/L 
     
  a is arbitrary integer value (preferably 2 to 4)
 
  Window (true space)   length= a L
  Window (adimensioned) length= a
  
  The equivalent cut-off length is 2*L.
  
  For interpolation :
      
    L=dX, and X = k dX, value at discrete nodes unchanged 
    
  
  For smoothing :
  
    L!=dX, otherwise value at discrete nodes unchanged (no smoothing)
    
  
  Weigths are compute at discrete node position:
  
  
  X= k dX, x=k dX/L
  
------------------------------------------------------------------------------*/
{
  double r,s,x;
  
  m=a*L/dx+1;
  
  weight=new T[m];
  
  for(int i=0;i<m;i++)
    weight[i]=0.;
 
  weight[0]=1;

  for(int i=1;i<m;i++) {
    x=M_PI*i*dx;
    r=sin(x)/x;
    x=M_PI*i*dx/a;
    s=sin(x)/x;
    weight[i]=r*s;
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Lanczos1D_BF(double L, double dx, int a, filter_t<float> & filter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=Lanczos1D_weights_template(L, dx, a, filter.weights, filter.nx);
  return(status);
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int Lanczos1D_BFatom_template(double L, double *x, int nvalues, int target, int a, T *buffer, T mask, T *BF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  
  Designed for series with irregular "space" sampling.
  
  No gain to split between weights computation and filtering (reversely to 
  regularly sampled sries).
  
------------------------------------------------------------------------------*/
{
  double r,s,xx,dx;
  int minf,msup,m;
  T w, sum;
  
  double x0=x[target];
  
  m=target;
  while(x[target]-x[m]<a*L) {
    if(m==0) break;
    m--;
    }
  minf=m;
  
  m=target;
  while(x[m]-x[target]<a*L) {
    if(m==nvalues-1) break;
    m++;
    }
  msup=m;
  
  m=MAX(target-minf,target-msup);
    
  BF[target]=buffer[target];
  sum=1.0;

  for(int i=minf;i<target;i++) {
    dx=(x[target]-x[i])/L;
    xx=M_PI*dx;
    r=sin(xx)/xx;
    xx=M_PI*dx/a;
    s=sin(xx)/xx;
    w=r*s;
    BF[target]+=w*buffer[i];
    sum+=w;
    }
      
  for(int i=target+1;i<msup+1;i++) {
    dx=(x[i]-x[target])/L;
    xx=M_PI*dx;
    r=sin(xx)/xx;
    xx=M_PI*dx/a;
    s=sin(xx)/xx;
    w=r*s;
    BF[target]+=w*buffer[i];
    sum+=w;
    }
  
  BF[target]/=sum;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int Lanczos1D_BF_template(double L, int a, double *x, int nvalues, T *buffer, T mask, T * & BF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  filter_t<float> filter;
  
  if(BF!=0) delete[] BF;
  
  BF=new T[nvalues];
  
  for(int n=0; n< nvalues; n++) {
    status=Lanczos1D_BFatom_template(L, x, nvalues, n, a, buffer, mask, BF);
    }
    
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Lanczos1D_BF(double L, int a, double *x, int nvalues, complex<float> *buffer, complex<float> mask, complex<float>* & BF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  filter_t<float> filter;
  
  status=Lanczos1D_BF_template(L, a, x, nvalues, buffer, mask, BF);
    
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int Lanczos1D_BF_template(double *L, int a, double *x, int nvalues, T *buffer, T mask, T * & BF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  filter_t<float> filter;
  
  if(BF!=0) delete[] BF;
  
  BF=new T[nvalues];
  
  for(int n=0; n< nvalues; n++) {
    status=Lanczos1D_BFatom_template(L[n], x, nvalues, n, a, buffer, mask, BF);
    }
    
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Lanczos1D_BF(double *L, int a, double *x, int nvalues, complex<float> *buffer, complex<float> mask, complex<float>* & BF)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  filter_t<float> filter;
  
  status=Lanczos1D_BF_template(L, a, x, nvalues, buffer, mask, BF);
    
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int Lanczos2D_weights(double L, double dx, double dy, int a, T * & weight, int & mx, int & my)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double r,s,x,xx,yy;
  
  mx=a*L/dx+1;
  my=a*L/dy+1;
  
  weight=new T[mx*my];
  
  for(int i=0;i<mx*my;i++)
    weight[i]=0.;
 
  weight[0]=1;

  for(int j=0;j<my;j++) {
    yy=j*dy;
    for(int i=0;i<mx;i++) {
      if((j==0) && (i==0)) continue;
      xx=i*dx;
      x=M_PI*sqrt(xx*xx+yy*yy);
      r=sin(x)/x;
      x=M_PI*i*dx/a;
      s=sin(x)/x;
      weight[mx*j+i]=r*s;
      }
    }  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Lanczos2D_weights(double L, double dx, double dy, int a, filter_t<float> & filter)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=Lanczos2D_weights(L, dx, dy, a,  filter.weights, filter.nx, filter.ny);
  return(status);
  
}

