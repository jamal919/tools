
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Damien Allain      LEGOS/CNRS, Toulouse, France

\brief tests equivalence between fftpack and gsl

Just for testing. The code is interesting to look at.
*/
/*----------------------------------------------------------------------------*/

#include "../config.h"

#if HAVE_LIBGSL == 1
#include <gsl/gsl_fft_real.h>
#endif

#include <stdio.h>

#if HAVE_LIBFFTPACK == 1
//dfft sont des routines de fftpack ....
extern "C" {
extern void dffti_(int *, double *);
extern void dfftf_(int *, double *,double *);
}
#endif

int main(void)
{
  int n,i,j;
  double *tf,*data,value;
  
  for(n=4;n<6;n++){
#if HAVE_LIBFFTPACK == 1
    tf=new double[n];
#endif
#if HAVE_LIBGSL == 1
    data=new double[n];
#endif
    
    for(j=0;j<n;j++){
      printf("vvv --- n=%d j=%d --- vvv\n",n,j);
      
      for(i=0;i<n;i++){
        value=(i==j);
        printf("%g\n",value);
        
#if HAVE_LIBFFTPACK == 1
        tf[i]=value;
#endif
#if HAVE_LIBGSL == 1
        data[i]=value;
#endif
        }
      printf("\n");
      
#if HAVE_LIBFFTPACK == 1
      {
      double *wsave=new double[2*n+15];
      dffti_(&n,         wsave);
      dfftf_(&n, &tf[0], wsave);
      delete[]wsave;
      }
#endif
      
#if HAVE_LIBGSL == 1
      {
      //see info:/gsl-ref/Mixed-radix FFT routines for real data
      gsl_fft_real_workspace *work = gsl_fft_real_workspace_alloc (n);
      gsl_fft_real_wavetable *wata = gsl_fft_real_wavetable_alloc (n);
      
      gsl_fft_real_transform (data, 1, n, wata, work);
      
      gsl_fft_real_wavetable_free (wata);
      gsl_fft_real_workspace_free (work);
      }
#endif

      for(i=0;i<n;i++){
#if HAVE_LIBFFTPACK == 1
        printf("%g ",tf[i]);
#endif
#if HAVE_LIBGSL == 1
        printf("%g ",data[i]);
#endif
        printf("\n");
        }
      printf("^^^ --- n=%d j=%d --- ^^^\n\n",n,j);
      }
    
#if HAVE_LIBFFTPACK == 1
    delete[]tf;
#endif
#if HAVE_LIBGSL == 1
    delete[]data;
#endif
    }
}