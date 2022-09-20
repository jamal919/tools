
/**************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Cyril Nguyen       LA/CNRS,    Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief FE integration optimisation
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h>
#include <string.h>
#include <cmath>

//#define MAIN_FE_INTEGRATION_SOURCE

#include "fe.h"
#include "geo.h"

int fe_chk_integrales(FILE *);

extern double (* fe_sproduct_ptrLGP1xLGP1_2D_a)  (double *, double *);
extern double (* fe_sproduct_ptrLGP1xLGP1_2D_b)  (double *, int);

extern double (* fe_sproduct_ptrLGP1xLGP1xLGP1_2D_a) (double *, double *, double *);
extern double (* fe_sproduct_ptrLGP1xLGP1xLGP1_2D_b) (double *, double *, int);

extern double (* fe_sproduct_ptrLGP1xLGP1xLGP1xLGP1_2D) (double *, double *, double *, int);

extern double (* fe_sproduct_ptrLGP1xNCP1_2D)      (double *, double *);
 
extern double (* fe_sproduct_ptrLGP1xLGP1xNCP1_2D_a)  (double *, double *, double *);
extern double (* fe_sproduct_ptrLGP1xLGP1xNCP1_2D_b)  (double *, int     , double *);
//
//
// /**-----------------------------------------------------------------------
// 2D nodal coordinates*/
// extern double fe_LGP1_x[3];
// extern double fe_LGP1_y[3];
//
// extern double fe_NCP1_x[3];
// extern double fe_NCP1_y[3];
//
// extern double fe_LGP2_x[6];
// extern double fe_LGP2_y[6];
//
//
// extern const double fe_LGP1xLGP1xNCP1_2D_c1;
// extern const double fe_LGP1xLGP1xNCP1_2D_c2;
// extern const double fe_LGP1xLGP1xNCP1_2D_c3;
// extern const double fe_LGP1xLGP1xNCP1_2D_c4;
//
// extern polynomial_t baseLGP2_x_2D[6];
// extern polynomial_t baseLGP2_y_2D[6];


//#define NOMINAL
#define OPTIMAL

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_optimize_integrales_R3(FILE *out,const char *name,double integrale[3][3][3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  int v,nvalue,npos[100],pos[100][100];
  double value[100];
  double mean=0,rms=0,count=0,error,z;

/* *----------------------------------------------------------------------------
  identify the identical integrale value*/
  nvalue=0;
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        m=9*i+3*j+k;
        for(v = 0; v < nvalue; v++) {
          if(fabs(value[v]-integrale[i][j][k]) < 1.e-06) {
            npos[v]++;
            pos[v][npos[v]-1]=m;
            goto next;
            }
          }
        nvalue++;
        value[nvalue-1]=integrale[i][j][k];
        npos[nvalue-1]=1;
        pos[nvalue-1][npos[nvalue-1]-1]=m;
next:
        continue;
        }
      }
    }

/* *----------------------------------------------------------------------------
  create C/C++ code for basic integrale computation*/
  for(v = 0; v < nvalue; v++) {
    n=0;
    m=pos[v][n];
    k=m%3;
    m=(m-k)/3;
    j=m % 3;
    m=(m-j)/3;
    i=m % 3;
//    fprintf(out,"\n tmp%d=fe_LGP1xLGP1xNCP1_2D[%d][%d][%d]*(p[%d]*q[%d]*r[%d]",nvalue,i,j,k,i,j,k);
    fprintf(out,"\n sum+=%s[%d][%d][%d]*(p[%d]*q[%d]*r[%d]",name,i,j,k,i,j,k);
    for(n = 1; n < npos[v]; n++) {
      m=pos[v][n];
      k=m%3;
      m=(m-k)/3;
      j=m % 3;
      m=(m-j)/3;
      i=m % 3;
      fprintf(out,"+p[%d]*q[%d]*r[%d]",i,j,k);
      }
    fprintf(out,") // %lf",value[v]);
    }
  fprintf(out,"\n");

/* *----------------------------------------------------------------------------
  create C/C++ code for integrale computation with last polynom being a base function*/
  for(k = 0; k < 3; k++) {
    nvalue=0;
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        m=9*i+3*j+k;
        for(v = 0; v < nvalue; v++) {
          if(fabs(value[v]-integrale[i][j][k]) < 1.e-06) {
            npos[v]++;
            pos[v][npos[v]-1]=m;
            goto next2;
            }
          }
        nvalue++;
        value[nvalue-1]=integrale[i][j][k];
        npos[nvalue-1]=1;
        pos[nvalue-1][npos[nvalue-1]-1]=m;
next2:
        continue;
        }
      }

    fprintf(out,"    case %d: \n",k);

    for(v = 0; v < nvalue; v++) {
      n=0;
      m=pos[v][n];
      k=m%3;
      m=(m-k)/3;
      j=m % 3;
      m=(m-j)/3;
      i=m % 3;
      fprintf(out,"\n      sum2+=%s[%d][%d][%d]*(p[%d]*q[%d]",name,i,j,k,i,j);
      for(n = 1; n < npos[v]; n++) {
        m=pos[v][n];
        k=m%3;
        m=(m-k)/3;
        j=m % 3;
        m=(m-j)/3;
        i=m % 3;
        fprintf(out,"+p[%d]*q[%d]",i,j);
//         if(n%3==0) {
//           fprintf(out,"\n");
//           for (i=0;i<strlen(name)+7+9+6;i++) fprintf(out," ");
//           }
        }
      fprintf(out,"); // %lf",value[v]);
      }
    fprintf(out,"\n");
    fprintf(out,"      break;\n\n");
    }

/* *----------------------------------------------------------------------------
  same, but slightly more optimal (no array indexing)*/
  for(k = 0; k < 3; k++) {
    nvalue=0;
    for(i = 0; i < 3; i++) {
      for(j = 0; j < 3; j++) {
        m=9*i+3*j+k;
        for(v = 0; v < nvalue; v++) {
          if(fabs(value[v]-integrale[i][j][k]) < 1.e-06) {
            npos[v]++;
            pos[v][npos[v]-1]=m;
            goto next3;
            }
          }
        nvalue++;
        value[nvalue-1]=integrale[i][j][k];
        npos[nvalue-1]=1;
        pos[nvalue-1][npos[nvalue-1]-1]=m;
next3:
        continue;
        }
      }

    fprintf(out,"    case %d: \n",k);

    for(v = 0; v < nvalue; v++) {
      n=0;
      m=pos[v][n];
      k=m%3;
      m=(m-k)/3;
      j=m % 3;
      m=(m-j)/3;
      i=m % 3;
      fprintf(out,"\n      sum2+=%.13e*(p[%d]*q[%d]",value[v],i,j);
      for(n = 1; n < npos[v]; n++) {
        m=pos[v][n];
        k=m%3;
        m=(m-k)/3;
        j=m % 3;
        m=(m-j)/3;
        i=m % 3;
        fprintf(out,"+p[%d]*q[%d]",i,j);
        }
      fprintf(out,"); // %lf",value[v]);
      }
    fprintf(out,"\n");
    fprintf(out,"      break;\n\n");
    }

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_optimize_integrales_R4(FILE *out, const char *name,double integrale[3][3][3][3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  int v,nvalue,npos[100],pos[100][100];
  double value[100];
  double mean=0,rms=0,count=0,error,z;

  nvalue=0;
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        for(l = 0; l < 3; l++) {
          m=27*i+9*j+3*k+l;
          for(v = 0; v < nvalue; v++) {
            if(fabs(value[v]-integrale[i][j][k][l]) < 1.e-06) {
              npos[v]++;
              pos[v][npos[v]-1]=m;
              goto next;
              }
            }
          nvalue++;
          value[nvalue-1]=integrale[i][j][k][l];
          npos[nvalue-1]=1;
          pos[nvalue-1][npos[nvalue-1]-1]=m;
next:
          continue;
          }
        }
      }
    }

  for(v = 0; v < nvalue; v++) {
    n=0;
    m=pos[v][n];
    l=m%3;
    m=(m-l)/3;
    k=m%3;
    m=(m-k)/3;
    j=m % 3;
    m=(m-j)/3;
    i=m % 3;
//    fprintf(out,"\n tmp%d=fe_LGP1xLGP1xNCP1_2D[%d][%d][%d]*(p[%d]*q[%d]*r[%d]",nvalue,i,j,k,i,j,k);
    fprintf(out,"\n sum+=%s[%d][%d][%d][%d]*(p[%d]*q[%d]*r[%d]*s[%d]",name,i,j,k,l,i,j,k,l);
    for(n = 1; n < npos[v]; n++) {
      m=pos[v][n];
      l=m%3;
      m=(m-l)/3;
      k=m%3;
      m=(m-k)/3;
      j=m % 3;
      m=(m-j)/3;
      i=m % 3;
      fprintf(out,"+p[%d]*q[%d]*r[%d]*s[%d]",i,j,k,l);
      if(n%3==0) {
        fprintf(out,"\n");
        for (i=0;i<strlen(name)+7+12;i++) fprintf(out," ");
        }
      }
    fprintf(out,"); // %lf",value[v]);
    }

  fprintf(out,"\n");

  for(l = 0; l < 3; l++) {
  nvalue=0;
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
          m=27*i+9*j+3*k+l;
          for(v = 0; v < nvalue; v++) {
            if(fabs(value[v]-integrale[i][j][k][l]) < 1.e-06) {
              npos[v]++;
              pos[v][npos[v]-1]=m;
              goto next2;
              }
            }
          nvalue++;
          value[nvalue-1]=integrale[i][j][k][l];
          npos[nvalue-1]=1;
          pos[nvalue-1][npos[nvalue-1]-1]=m;
next2:
          continue;
          }
        }
      }

    fprintf(out,"    case %d: \n",l);

    for(v = 0; v < nvalue; v++) {
      n=0;
      m=pos[v][n];
      l=m%3;
      m=(m-l)/3;
      k=m%3;
      m=(m-k)/3;
      j=m % 3;
      m=(m-j)/3;
      i=m % 3;
//    fprintf(out,"\n tmp%d=fe_LGP1xLGP1xNCP1_2D[%d][%d][%d]*(p[%d]*q[%d]*r[%d]",nvalue,i,j,k,i,j,k);
      fprintf(out,"\n      sum+=%s[%d][%d][%d][%d]*(p[%d]*q[%d]*r[%d]",name,i,j,k,l,i,j,k);
      for(n = 1; n < npos[v]; n++) {
        m=pos[v][n];
        l=m%3;
        m=(m-l)/3;
        k=m%3;
        m=(m-k)/3;
        j=m % 3;
        m=(m-j)/3;
        i=m % 3;
        fprintf(out,"+p[%d]*q[%d]*r[%d]",i,j,k);
        if(n%3==0) {
          fprintf(out,"\n");
          for (i=0;i<strlen(name)+7+12+6;i++) fprintf(out," ");
          }
        }
      fprintf(out,"); // %lf",value[v]);
      }
    fprintf(out,"\n");
    fprintf(out,"      break;\n\n");

    }
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_chk_integrales(FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n,status;
  int v,nvalue,npos[100],pos[100][100];
  double value[100];
  double p[6], q[6];
  double LG[3][3][2];
  double x,y,nx,ny,L[3],angle;
  double sum;
  double beta[3][3];
  double beta_x[3], beta_y[3];
  double beta_t[3], beta_p[3];
  double gamma_x[3], gamma_y[3];
  double gamma_t[3], gamma_p[3];
  double mean=0,rms=0,count=0,error,z;
  double z0,z1,z2,z3;

  for(i=0;i<3;i++) {
    j=i;
    LG[i][j][0]=(double)  1.;
    LG[i][j][1]=(double)  1.;
    j=(i+1) % 3;
    LG[i][j][0]=(double)  1.;
    LG[i][j][1]=(double) -1.;
    j=(i+2) % 3;
    LG[i][j][0]=(double) -1.;
    LG[i][j][1]=(double)  1.;
    }
/*-----------------------------------------------------------------------
  LGP1 base function on NCP1 node*/
  x = 0.5;
  y = 0.5;
  fe_LGP1base(x, y, beta[0]);
  x = 0.0;
  y = 0.5;
  fe_LGP1base(x, y, beta[1]);
  x = 0.5;
  y = 0.0;
  fe_LGP1base(x, y, beta[2]);

/*-----------------------------------------------------------------------
  derivative of LGP1 base function on reference triangle*/
  x = 0.;
  y = 0.;
  fe_LGP1prime(x, y, beta_x,  beta_y);

/*-----------------------------------------------------------------------
  derivative of NCP1 base function on reference triangle*/
  fe_NCP1prime(x, y, gamma_x, gamma_y);

  for(i=0;i<3;i++) {
    p[i]=1.0;
    q[i]=1.0;
    }
  sum=fe_integraleLGP1xLGP1_2D(p,q);
  for(i=0;i<6;i++) {
    p[i]=1.0;
    }
  sum=fe_integraleLGP2_2D(p);

  for(i=0;i<3;i++) {
    sum=0;
    for(j=0;j<3;j++) {
      if(i==j) {
        sum+=p[i];
        }
      }
    }

  nvalue=0;
  for(i = 0; i < 3; i++) {
    for(j = 0; j < 3; j++) {
      for(k = 0; k < 3; k++) {
        m=9*i+3*j+k;
        for(v = 0; v < nvalue; v++) {
          if(fabs(value[v]-fe_LGP1xLGP1xNCP1_2D[i][j][k]) < 1.e-06) {
            npos[v]++;
            pos[v][npos[v]-1]=m;
            goto next;
            }
          }
        nvalue++;
        value[nvalue-1]=fe_LGP1xLGP1xNCP1_2D[i][j][k];
        npos[nvalue-1]=1;
        pos[nvalue-1][npos[nvalue-1]-1]=m;
next:
        continue;
        }
      }
    }

  for(v = 0; v < nvalue; v++) {
    n=0;
    m=pos[v][n];
    k=m%3;
    m=(m-k)/3;
    j=m % 3;
    m=(m-j)/3;
    i=m % 3;
//    fprintf(out,"\n tmp%d=fe_LGP1xLGP1xNCP1_2D[%d][%d][%d]*(p[%d]*q[%d]*r[%d]",nvalue,i,j,k,i,j,k);
    fprintf(out,"\n sum+=fe_LGP1xLGP1xNCP1_2D[%d][%d][%d]*(p[%d]*q[%d]*r[%d]",i,j,k,i,j,k);
    for(n = 1; n < npos[v]; n++) {
      m=pos[v][n];
      k=m%3;
      m=(m-k)/3;
      j=m % 3;
      m=(m-j)/3;
      i=m % 3;
      fprintf(out,"+p[%d]*q[%d]*r[%d]",i,j,k);
      }
    fprintf(out,") // %lf",value[v]);
    }
  fprintf(out,"\n");

  status=fe_optimize_integrales_R3(out,"fe_LGP1xLGP1xLGP1_2D",fe_LGP1xLGP1xLGP1_2D);
  status=fe_optimize_integrales_R3(out,"fe_LGP1xLGP1xNCP1_2D",fe_LGP1xLGP1xNCP1_2D);
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Note:

  19/10/2009

    serious bug found : wrong optimal formula taken from wrong
    fe_optimize_integrales_R3 call:

  status=fe_optimize_integrales_R3(out,"fe_LGP1xNCP1xNCP1_2D",fe_LGP1xLGP1xNCP1_2D);

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  status=fe_optimize_integrales_R3(out,"fe_LGP1xNCP1xNCP1_2D",fe_LGP1xNCP1xNCP1_2D);
  status=fe_optimize_integrales_R4(out,"fe_LGP1xLGP1xLGP1xLGP1_2D",fe_LGP1xLGP1xLGP1xLGP1_2D);
  status=fe_optimize_integrales_R4(out,"fe_LGP1xLGP1xNCP1xNCP1_2D",fe_LGP1xLGP1xNCP1xNCP1_2D);

  return (0);
}

