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

#include "tools-define.h"
#include "tools-structures.h"

#include "matrix.h"

#include "geo.h"
#include "map.h"
#include "grd.h"
#include "bmg.h"
#include "polygones.h"
#include "tides.h"
#include "tides.def"
#include "filter.h"
#include "statistic.h"
#include "mss.v2.h"
#include "rutin.h"

double DEG2RAD = M_PI/180.0;
double RAD2DEG = 180.0/M_PI;

double rTerre = 6378.1363e3;
double omegaTerre = 7.292115e-5;
double j2Terre = 1082.63622070e-6;
double fTerre = 298.257;
//double demiGdAxe = 7714.43e3;
double demiGdAxe = 7714.4278e3;
double moyenMouv = 2.0*M_PI/6745.72;
//double inclin = 66.04*DEG2RAD;
double inclin = 66.0075*DEG2RAD;
double omegaP = -3.0 / 2.0 * moyenMouv * j2Terre * ( rTerre / demiGdAxe )*( rTerre / demiGdAxe ) * cos( inclin );



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void trace( int nt, double *tmtna, double lambda0, double *lon, double *lat)
{
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  //   Routine de calcul de la trace au sol d'un satellite en orbite circulaire
  //
  // * Arguments
  //
  //   nt      (E) : nombre de dates
  //   tmtna   (E) : duree entre le point courant et le noeud ascendant (en s)
  //   lambda0 (E) : longitude du noeud ascendant (en rad)
  //   lon     (S) : latitudes de la trace (en rad)
  //   lat     (S) : longitudes de la trace (en rad)
  //

  double phi, lambdaX, lambdaY;
  int i;

  for(i=0;i<nt;i++) {

      phi= asin( sin( moyenMouv * tmtna[i] ) * sin( inclin ) ) ;
      lambdaX = acos( (cos(lambda0+(omegaP-omegaTerre)*tmtna[i])*cos(moyenMouv*tmtna[i]) - sin(lambda0+(omegaP-omegaTerre)*tmtna[i])*sin(moyenMouv*tmtna[i])*cos(inclin)) / cos(phi) );
      lambdaY = asin( (sin(lambda0+(omegaP-omegaTerre)*tmtna[i])*cos(moyenMouv*tmtna[i]) + cos(lambda0+(omegaP-omegaTerre)*tmtna[i])*sin(moyenMouv*tmtna[i])*cos(inclin)) / cos(phi) );
      //calcul de la longitude
      if( lambdaX <  M_PI/2.0 ){if( lambdaY <  0.0 )lon[i] = lambdaY;if( lambdaY >= 0.0 )lon[i] = lambdaX;}
      if( lambdaX >= M_PI/2.0 ){if( lambdaY <  0.0 )lon[i] = M_PI-lambdaY;if( lambdaY >= 0.0 )lon[i] = lambdaX;}
      //calcul de la latitude
      lat[i] = atan( tan(phi)/((1.0-1.0/fTerre)*(1.0-1.0/fTerre)) );

    } //fin boucle for i<nt

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void interpTrace( char tTrace[3], double lonRef, double latRef, int nt, double *dt, double *lon, double *lat )
{
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  //   Routine de calcul des points de la trace au sol, a partir d'un point de reference.
  //
  // * Argument
  //
  //   nt (E)     : nombre de dates
  //   dt (E)     : duree entre le point reference et les points interpoles
  //                de la trace au sol (en s)
  //   latRef (E) : latitude du point de reference (en rad)
  //   lonRef (E) : longitude du point de reference (en rad)
  //   tTrace (E) : type de trace ('ASC' ou 'DSC')
  //   lat (S)    : latitude des points interpoles de la trace (en rad)
  //   lon (S)    : longitude des points interpoles de la trace (en rad)
  int i;
  double  phiRef;
  double  xSat, ySat, ksi;
  double  coefA, coefB, cBscA;
  double  cpsi, spsi, psi;
  double  lambda0, corrLon, corrLat;

  double tmtna,lonCalc,latCalc;

  phiRef    = atan( tan( latRef ) * ( 1.0 - 1.0/fTerre)*( 1.0 - 1.0/fTerre) );

  //Calcul de tmtna
  tmtna  = asin( sin( phiRef )/sin( inclin ) ) / moyenMouv;    //temps ecoule sur la trace depuis le passage a l'equateur de la trace ascendante!!!
  if(strcmp(tTrace,"ASC")!=0)tmtna = M_PI/moyenMouv - tmtna;   //donc on part du passage a l'equateur de la trace descendante (M_PI/moyenMouv) et on revient a la bonne latitude!!!

  //Calcul de lambda0
  xSat  = cos( phiRef )*cos( lonRef );
  ySat  = cos( phiRef )*sin( lonRef );

  ksi   = moyenMouv*tmtna;    //ksi = angle effectue sur le cercle de l'orbite du satellite
  coefA = cos( ksi );
  coefB = sin( ksi )*cos( inclin );
  if(strcmp(tTrace,"ASC")!=0)coefB*=-1.0;   //car l'inclinaison apparente sur une trace descendante est -66.00 et pas 66.00

  cBscA = coefB/coefA;

  spsi  = (ySat-xSat*cBscA)/(coefA+cBscA);
  cpsi  = (xSat+coefB*spsi)/coefA;

  psi   = asin( spsi );
  if ( cpsi < 0.0 ) psi = M_PI - psi;
  lambda0 = psi - (omegaP-omegaTerre)*tmtna;
  
  //Calcul des corrections

  trace( 1, &tmtna, lambda0, &lonCalc, &latCalc );

  corrLon = lonCalc - lonRef;
  corrLat = latCalc - latRef;
  //Calcul des points interpoles

  for(i=0;i<nt;i++)dt[i]+=tmtna;
  trace( nt, dt, lambda0, lon, lat );

  for(i=0;i<nt;i++) {
      lat[i] = lat[i] - corrLat;
      lon[i] = lon[i] - corrLon;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int compute_nominal_track(grid_t *backbone,grid_t *nominal_track)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  //backbone contient la meilleure trace nominal possible
  //donc on va chercher dans backbone un point tous les 7 km

  //pour les projections on va utiliser lib_proj4

  double lon1,lon2,lat1,lat2,ref_lon,ref_lat,interval,distance,restant,longeur;
  int n,test,av;
  double x,y;
  //char *parms[] = { "proj=lcc", "ellps=GRS80", "lon_0=0E"};
  projPJ ref;
  projUV xy1,xy2,xy;
  projUV lp1,lp2,lp;

  const char *parms = "proj=merc +ellps=GRS80 +lon_0=0E";

  if ( ! (ref = pj_init_plus(parms)) ) {
    __ERR_BASE_LINE__( "Projection initialization failed\n"); exit(1);
    }

  interval=5.9;  //interval entre les points de la trace nominale

  longeur=geo_km(backbone->x[0],backbone->y[0],backbone->x[backbone->ny-1],backbone->y[backbone->ny-1]);
  longeur/=5.9;

  nominal_track->x=(double *)calloc((int)(longeur+100),sizeof(double));
  nominal_track->y=(double *)calloc((int)(longeur+100),sizeof(double));

  nominal_track->x[0]=backbone->x[0];
  nominal_track->y[0]=backbone->y[0];

  test=0;
  av=0;
  distance=0;
  n=0;
  while(test==0) {//cette boucle ne se termine jamais ...
    if(n==backbone->ny)goto cleanUp;
    while(distance<interval){
      restant=interval-distance;
      longeur=geo_km(backbone->x[n],backbone->y[n],backbone->x[n+1],backbone->y[n+1]);
      distance+=longeur;
      n++;
      if(n>backbone->ny)goto cleanUp;
      }

    //la solution se trouve entre n et n-1
    //donc on calcule la position en KM
    //il faut donc passer en projection ::

    lp1.u=backbone->x[n-1]*DEG_TO_RAD;
    lp1.v=backbone->y[n-1]*DEG_TO_RAD;
    lp2.u=backbone->x[n]*DEG_TO_RAD;
    lp2.v=backbone->y[n]*DEG_TO_RAD;

    xy1=pj_fwd(lp1, ref);
    xy2=pj_fwd(lp2, ref);

    xy.u=xy1.u+(xy2.u-xy1.u)*restant/longeur;
    xy.v=xy1.v+(xy2.v-xy1.v)*restant/longeur;

    lp=pj_inv(xy,ref);
    lp.u*=RAD_TO_DEG;
    lp.v*=RAD_TO_DEG;

    av++;
    nominal_track->x[av]=lp.u;
    nominal_track->y[av]=lp.v;
    nominal_track->nx=av;
    nominal_track->ny=av;

    distance=geo_km(lp.u,lp.v,backbone->x[n],backbone->y[n]);
    }

cleanUp:
  pj_free(ref);
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_topographic(double *Cm, int dimM, grid_t mssgrid, float *mss_cls)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*-----------------------------------------------------------------------------
  simple test : set gaussian covariances*/
  int i,j,k,l,m;
  int rstatus;
  float  *P2,*D,dx2;
  double x,y;
  float *vector,dum;
  float z,mask;
  float *topo,*topox,*topoy,*tmp;
  char input[1024],output[1024];
  grid_t topogrid3d;

  double norm;

/*------------------------------------------------------------------------------
  load bottom topography */

  sprintf(input,"/data/ocean/produits/topography/data/albicocca.grd");

  rstatus=grd_loadgrid(input,&topogrid3d);
  //rstatus=map2d_ completegridaxis(&topogrid);
  rstatus=map3d_completegridaxis(&topogrid3d);

  topo=(float *)malloc(topogrid3d.nx*topogrid3d.ny*sizeof(float));
  rstatus= grd_loadr1(input,topogrid3d,topo,&mask);

/*-----------------------------------------------------------------------------
  surface curvature constraint

  minimize a weighted (s")^2

  sxC s = (s")^2 = (sP")x D P"s

where:

  D  : diagonale weighting matrix
  P" : second derivative operator

  P" is a highly banded matrix

-----------------------------------------------------------------------------*/

  topox=(float *)malloc(topogrid3d.nx*topogrid3d.ny*sizeof(float));
  topoy=(float *)malloc(topogrid3d.nx*topogrid3d.ny*sizeof(float));

//  rstatus=copy_grid_to_grid1d(topogrid3d,&topogrid1d);

  rstatus=map_curve(topogrid3d, topo, mask,(int) 0, topox,topoy);

  sprintf(output,"../data/topo.crv");
  rstatus=bmg_saver2(output, 1, 1, 1, topogrid3d, topox, topoy, (float) 0., mask);


/*-----------------------------------------------------------------------------
  surface seconde derivative coefficients:
  ----------------------------------------

  1st order approximation:

  s"=[s(m+1)+s(m-1)-2*s(m)]/dx^2

  D is a modified second derivative
  equivalent to a correlation, or an inverse L^2/range

  range=0.1 metres, L=0.1deg D=10.


-----------------------------------------------------------------------------*/

  P2=(float *)malloc(dimM*dimM*sizeof(float));
  D=(float *)malloc(dimM*sizeof(float));

  for(k=0;k<dimM*dimM;k++) P2[k]=0;
  for(k=0;k<dimM;k++) D[k]=0;

  for(j=0;j<mssgrid.ny;j++)
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      x=mssgrid.x[k];
      y=mssgrid.y[k];
      rstatus=map_interpolation(topogrid3d, topox,mask,x,y,&z);
      if (z!=mask) D[k]=z*z;
      rstatus=map_interpolation(topogrid3d, topoy,mask,x,y,&z);
      if (z!=mask) D[k]+=z*z;
/*-----------------------------------------------------------------------------
      appropriate if unknown is mss itself, not mss error*/
      D[k]=10*sqrt(D[k])*mss_cls[k]/fabs (topo[k]);
      D[k]=1000.;
      }

/*-----------------------------------------------------------------------------
  compute D x P2 along x direction directly in P2*/

  dx2=mssgrid.dx*mssgrid.dx;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      P2[k*dimM+k]=-2./dx2/D[k];
      if(i< mssgrid.nx-1) P2[k*dimM+k+1]=1./dx2/D[k];
      if(i>0)             P2[k*dimM+k-1]=1./dx2/D[k];
      }
    }

/*-----------------------------------------------------------------------------
  compute P2x x D x P2  along x direction

   ( k, l) =   [D x P2]x (k,m)  x [D x P2]  (m,l)
   ( k, l) =   [D x P2]  (m,k)  x [D x P2]  (m,l)

  it is k col of [D x P2] x l col of [D x P2] product
-----------------------------------------------------------------------------*/

  printf("upper bw of P2 %d\n",gmatrix_ubw(P2,dimM,dimM));
  printf("lower bw of P2 %d\n",gmatrix_lbw(P2,dimM,dimM));

  vector =(float *)malloc(dimM*sizeof(float));
  for(l=0;l<dimM;l++) {
    for(m=0;m<dimM;m++) vector[m]=P2[l*dimM+m];

#if   ATLAS == 1
    norm=cblas_dsdot(dimM, vector, 1, vector, 1);
#elif CBLAS == 1
    norm=cblas_dsdot(dimM, vector, 1, vector, 1);
#else
#error "stuff missing"
#endif

    if(norm==0.) {
      continue;
      }
    for(k=MAX(0,l-4);k<=MIN(dimM-1,l+4);k++) {
#if ATLAS == 1
      dum=cblas_dsdot(dimM, &(P2[l*dimM]), 1, &(P2[k*dimM]), 1);
#elif CBLAS == 1
      dum=cblas_dsdot(dimM, &(P2[l*dimM]), 1, &(P2[k*dimM]), 1);
#else
#error "stuff missing"
#endif
      Cm[l*dimM+k]+=dum;
      }
    }
/*-----------------------------------------------------------------------------
  internal consistency chek, Cm must be symmetric */
  rstatus=dgmatrix_symmetry(Cm, dimM);

/*-----------------------------------------------------------------------------
  compute D x P2 along y direction directly in P2*/

  for(k=0;k<dimM*dimM;k++) P2[k]=0;

  dx2=mssgrid.dy*mssgrid.dy;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      P2[k*dimM+k]=-2./dx2/D[k];
      if(j< mssgrid.ny-1) P2[k*dimM+k+dimM]=1./dx2/D[k];
      if(j>0)             P2[k*dimM+k-dimM]=1./dx2/D[k];
      }
    }

/*-----------------------------------------------------------------------------
  compute P2x x D x P2  along y direction

   ( k, l) =   [D x P2]x (k,m)  x [D x P2]  (m,l)
   ( k, l) =   [D x P2]  (m,k)  x [D x P2]  (m,l)

  it is k col of [D x P2] x l col of [D x P2] product
-----------------------------------------------------------------------------*/

  for(l=0;l<dimM;l++) {
    for(m=0;m<dimM;m++) vector[m]=P2[l*dimM+m];
#if ATLAS == 1
    norm=cblas_dsdot(dimM, vector, 1, vector, 1);
#elif CBLAS == 1
    norm=cblas_dsdot(dimM, vector, 1, vector, 1);
#else
#error "stuff missing"
#endif
    if(norm==0.) {
      continue;
      }
    for(k=MAX(0,l-4);k<=MIN(dimM-1,l+4);k++) {
#if ATLAS == 1
      dum=cblas_dsdot(dimM, &(P2[l*dimM]), 1, &(P2[k*dimM]), 1);
#elif CBLAS == 1
      dum=cblas_dsdot(dimM, &(P2[l*dimM]), 1, &(P2[k*dimM]), 1);
#else
#error "stuff missing"
#endif
      Cm[l*dimM+k]+=dum;
      }
    }

  free(vector);vector=NULL;
  free(tmp);tmp=NULL;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_gaussian(double *Cm, int dimM, grid_t mssgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*-----------------------------------------------------------------------------
  simple test : set gaussian covariances*/

  double *A,*B,*C1,*C2;
  int     neq,job,*pivot,rstatus,nrhs;
  int     i,j,k,l,m,n,kk,ll;
  double *window,lengthx,lengthy;
  double  dx,dy;
  int     hlwx,hlwy,lwx,lwy;


  lwx=11;
  lwy=21;

  hlwx=(lwx-1)/2;
  hlwy=(lwy-1)/2;

  window=(double *)malloc(lwx*lwy*sizeof(double));

/*-----------------------------------------------------------------------------
  simple test : set gaussian covariances*/
  lengthx=0.025;
  lengthy=0.025;

/*   lengthx=0.02; */
/*   lengthy=0.1; */

  for(j=0;j<lwy;j++) {
    for(i=0;i<lwx;i++) {
      k=j*lwx+i;
      dx=(i-hlwx)*mssgrid.dx/lengthx;
      dy=(j-hlwy)*mssgrid.dy/lengthy;
      window[k]=0.1*exp(-0.5*(dx*dx+dy*dy));
/*       window[k]=0.01*1./(1+(dx*dx+dy*dy)); */
      }
    }

  A=(double *)malloc(dimM*dimM*sizeof(double));
  B=(double *)malloc(dimM*dimM*sizeof(double));
  C2=(double *)malloc(dimM*dimM*sizeof(double));
  for(k=0;k<dimM*dimM;k++) A[k]=0;
  for(k=0;k<dimM*dimM;k++) C2[k]=0;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      for(ll=MAX(j-hlwy,0);ll<=MIN(mssgrid.ny-1,j+hlwy);ll++) {
        for(kk=MAX(i-hlwx,0);kk<=MIN(mssgrid.nx-1,i+hlwx);kk++) {
          l=ll*mssgrid.nx+kk;
          A[l*dimM+k]=window[(ll-j+hlwy)*lwx+(kk-i+hlwx)];
          }
        }
      }
    }

/*-----------------------------------------------------------------------------
  internal consistency check, A must be symmetric */
  rstatus=dgmatrix_symmetry(A, dimM);

  rstatus= dgmatrix_transpose(A, B, dimM, dimM);
  C1=dmatrix_product(A, B, dimM, dimM, dimM);

  neq=dimM;

  goto skip_test;

/*-----------------------------------------------------------------------------
  internal consistency check, C1 must be symmetric and positive (+definite ?) */


#if ATLAS == 1
    rstatus=clapack_dpotrf(CblasRowMajor, CblasUpper, neq,C1,neq);
#elif LAPACKF_ == 1
//    rstatus=dpotrf_(CblasRowMajor, CblasUpper, &neq, C1, &neq);
#else
#error "stuff missing"
#endif
  if(rstatus==0) printf("C1 found to be symmetric and positive; good job!!!\n");
  else {
    __OUT_BASE_LINE__("C1 found to be non-symmetric or non-positive (rstatus=%d); bad boy!!!\n",rstatus);
    exit(-1);
    }

/*-----------------------------------------------------------------------------
  C1 has been factored to perform tests, we must rebuild it*/
    free(C1);C1=NULL;
  C1=dmatrix_product(A, B, dimM, dimM, dimM);

 skip_test:

/*-----------------------------------------------------------------------------
  internal consistency chek, C1 must be symmetric */
  rstatus=dgmatrix_symmetry(C1, dimM);

/*-----------------------------------------------------------------------------
  inverse the covariances*/
  pivot =(int *)malloc(dimM*sizeof(int));
  job=0;

  rstatus=poc_getrf(neq, C1, pivot);
  for(l=0;l<dimM;l++) {
     /*     for(k=0;k<dimM;k++) C2[l*dimM+k]=0.; */
    C2[l*dimM+l]=1.;
    rstatus=poc_getrs(neq,nrhs,C1,pivot,&C2[l*dimM]);
    }

/*-----------------------------------------------------------------------------
  internal consistency chek, Cm must be symmetric */
  rstatus=dgmatrix_symmetry(C2, dimM);

/*-----------------------------------------------------------------------------
  copy C2 in Cm, in such a way that Cm is built to be symmetric */
  for (l=0;l<dimM;l++)
    for (k=l;k<dimM;k++) {
    m=dimM*k+l;
    n=dimM*l+k;
    Cm[n]=0.5*(C2[m]+C2[n]);
    Cm[m]=Cm[n];
    }

  free(C1);C1=NULL;
  free(A);A=NULL;
  free(B);B=NULL;
  free(C2);C2=NULL;
  free(pivot);pivot=NULL;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int set_byderivative(double *Cm, int dimM, grid_t mssgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*-----------------------------------------------------------------------------
  s*/

  double *A,*B,*C1,*C2;
  int     neq,job,*pivot,rstatus;
  int     i,j,k,l,m,n,kk,ll;
  double *window,lengthx,lengthy;
  double  dx,dy,x,y,variance;
  double  dx2,dy2;

/*-----------------------------------------------------------------------------
  correlation distance*/
  lengthx=0.1;
  lengthy=0.05;


/*-----------------------------------------------------------------------------
  normalize*/
  lengthx/=mssgrid.dx;
  lengthy/=mssgrid.dy;

  dx2=lengthx*lengthx;
  dy2=lengthy*lengthy;

  variance=0.1;

  A=(double *)malloc(dimM*dimM*sizeof(double));
  B=(double *)malloc(dimM*dimM*sizeof(double));
  
  for(k=0;k<dimM*dimM;k++) A[k]=0;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      x=mssgrid.x[k];
      y=mssgrid.x[k];

      if((i!=0) && (i!=mssgrid.nx-1)) {
        l=j*mssgrid.nx+i-1;
        A[l*dimM+k]=dx2;
        l=j*mssgrid.nx+i+1;
        A[l*dimM+k]=dx2;
        A[k*dimM+k]-=2.*dx2;
        }

/*-----------------------------------------------------------------------------
      shifted condition*/
/*
      if(i==0) {
        l=j*mssgrid.nx+i;
        A[l*dimM+k]=dx2*100;
        l=j*mssgrid.nx+i+2;
        A[l*dimM+k]=dx2*100;
        l=j*mssgrid.nx+i+1;
        A[k*dimM+k]-=2*dx2*100;
        }
      if(i==mssgrid.nx-1) {
        l=j*mssgrid.nx+i-2;
        A[l*dimM+k]=dx2*100;
        l=j*mssgrid.nx+i;
        A[l*dimM+k]=dx2*100;
        l=j*mssgrid.nx+i-1;
        A[k*dimM+k]-=2*dx2*100;
        }

*/
/*-----------------------------------------------------------------------------
      continuty condition*/
/*

      if(i==0) {
        l=j*mssgrid.nx+i+1;
        A[l*dimM+k]=dx2;
        A[k*dimM+k]-=2*dx2;
        }
      if(i==mssgrid.nx-1) {
        l=j*mssgrid.nx+i-1;
        A[l*dimM+k]=dx2;
        A[k*dimM+k]-=2*dx2;
        }
*/
/*-----------------------------------------------------------------------------
  slope condition*/

      if(i==0) {
        l=j*mssgrid.nx+i;
        A[l*dimM+k]=-1./mssgrid.dx;
        l=j*mssgrid.nx+i+1;
        A[l*dimM+k]=+1./mssgrid.dx;
        }

      if(i==mssgrid.nx-1) {
        l=j*mssgrid.nx+i-1;
        A[l*dimM+k]=+1./mssgrid.dx;
        l=j*mssgrid.nx+i;
        A[l*dimM+k]=-1./mssgrid.dx;
        }


      }
    }

  rstatus= dgmatrix_transpose(A, B, dimM, dimM);
  C1=dmatrix_product02(A, B, dimM, dimM, dimM);

  for(k=0;k<dimM*dimM;k++) A[k]=0;

  for(j=0;j<mssgrid.ny;j++) {
    for(i=0;i<mssgrid.nx;i++) {
      k=j*mssgrid.nx+i;
      x=mssgrid.x[k];
      y=mssgrid.x[k];

      if((j!=0) && (j!=mssgrid.ny-1)) {
        l=(j-1)*mssgrid.nx+i;
        A[l*dimM+k]=dy2;
        l=(j+1)*mssgrid.nx+i;
        A[l*dimM+k]=dy2;
        A[k*dimM+k]-=2*dy2;
        }

      }
    }

  rstatus= dgmatrix_transpose(A, B, dimM, dimM);
  C2=dmatrix_product02(A, B, dimM, dimM, dimM);

/*-----------------------------------------------------------------------------
  internal consistency check, A must be symmetric */
  rstatus=dgmatrix_symmetry(C1, dimM);
  rstatus=dgmatrix_symmetry(C2, dimM);

/*-----------------------------------------------------------------------------
  copy A in Cm, in such a way that Cm is built to be symmetric */

  for(k=0;k<dimM*dimM;k++) Cm[k]=0;

  for (l=0;l<dimM;l++) {
    m=dimM*l+l;
    Cm[m]=1./variance;
    }

  for (l=0;l<dimM;l++)
    for (k=l;k<dimM;k++) {
    m=dimM*l+k;
    Cm[m]+=C1[m]+C2[m];
    }

  for (l=0;l<dimM;l++)
    for (k=l;k<dimM;k++) {
    m=dimM*k+l;
    n=dimM*l+k;
    Cm[m]=Cm[n];
    }

  free(A);A=NULL;
  free(B);B=NULL;
  free(C1);C1=NULL;
  free(C2);C2=NULL;

  return(0);
}




//#####################################################################################################
//
//fonction d'inversion des donnes etc etc ...
//
//
//#####################################################################################################

int compute(grid_t mssgrid,serie_t TPdata,float *mss_cls,double *error,float *innovation,float *analysis,float *G)
{
  double  *A,*b;
  int dimM,dimD,m,n,l,count;
  double *Cm,weights[4];
  int k,ldaD,diagD;
  float  *Cd,*tG,*tmp1,*tmp2,*tmp3,*model,*observation;
  double x,y;
  int rstatus;
  int    nodes[4];
  float dum;
  float  *variance1;
  int nindex,*index,neq,job,*pivot;
  float  *svector,snorm;
  double *dvector,dnorm;


/*-----------------------------------------------------------------------------

Inverse problem solution
------------------------

For details and foramlism, see Tarantola (Elsevier), pp197

Penalty function J:

J=(d-Gm)xCddeg(d-Gm)+(m-m0)xCmdeg(m-m0)

Minimizing J (analysis step):

m(optimal)=m'=m0+[GxxCddegxG+Cmdeg]degxGxxCddeg(d-Gm)

Posterior covariance of inverse analysis:

Cm'=[GxxCddegxG+Cmdeg]deg

where:

deg denotes inverse
x denotes transpose

Cd  = prior error covariance on data
Cm  = prior error covariance on model
Cm' = posterior error covariance on analysis
G   = interpolation/observation operator

Cddeg=inverse of Cd
Cmdeg=inverse of Cm
Gx = transpose of G

Cd is dimDxdimD square matrix (huge in this application)
Cm is dimMxdimM square matrix (reasonable in this application)

GxxCddeg is dimMxdimG rectangular matrix

I. Quick analysis problem
xxxxxxxxxxxxxxxxxxxxxxxx

3 steps technique:

1) form GxxCddeg and [GxxCddegxG+Cmdeg]

2) solve [GxxCddegxG+Cmdeg]X=GxxCddeg(d-Gm)

3) m'=m0+X

posterior covariances not available !!!
---------------------------------------

II. Full analysis problem
xxxxxxxxxxxxxxxxxxxxxxxx

3 steps technique:

1) form GxxCddeg and [GxxCddegxG+Cmdeg]

2) explicitely compute Cm'=[GxxCddegxG+Cmdeg]deg

3) m'=m0+Cm'x GxxCddeg(d-Gm)

posterior covariances available !!!
-----------------------------------

III. Consistency check
xxxxxxxxxxxxxxxxxxxxxx

Jminimum=(d-Gm')xCddeg(d-Gm')+(m'-m0)xCmdeg(m'-m0) ~ Nobs

-----------------------------------------------------------------------------*/
/*   TPdata.count=1; */
/*-----------------------------------------------------------------------------
  matrix leading dimensions */
  dimM=mssgrid.nx*mssgrid.ny;
  dimD=TPdata.count;

/*-----------------------------------------------------------------------------
  square inverse covariance matrix initialisation.
  Cm and Cd are actually the inverse matrices */
  Cm=(double *)malloc(dimM*dimM*sizeof(double));
  for(k=0;k<dimM*dimM;k++) Cm[k]=0.;

/*-----------------------------------------------------------------------------
  Cd is (horizontal) band matrix; diagonal at row 1+ldaD/2 */
  ldaD=1;
  diagD=ldaD/2;
  printf("peut etre un BUG\n");
  Cd=(float *)malloc(dimD*ldaD*sizeof(float));
  for(k=0;k<dimD*ldaD;k++) Cd[k]=0.;

/*-----------------------------------------------------------------------------
  diagonal (inverse variance) constraint */
  for(k=0;k<dimM;k++) Cm[k*dimM+k]=1./0.01; /* 100 cm^2 =0.01 m^2 */
  for(k=0;k<dimD;k++) Cd[k*ldaD+diagD]=1./0.002; /* 20 cm^2 =0.002 m^2 */

/*-----------------------------------------------------------------------------
  observation matrix initialisation */
  if( (tG=(float *)malloc(dimM*dimD*sizeof(float))) == NULL) gmerror("can not allocate tG matrix");
  for(k=0;k<dimM*dimD;k++) G[k]=0;

/*-----------------------------------------------------------------------------
  get data point interpolation coeffcients (observation operator)*/
  for(k=0;k<dimD;k++) {
    x=TPdata.data[k].lon;
    y=TPdata.data[k].lat;
    rstatus=map_coeffcients02(mssgrid,x,y,weights,nodes,&count);
    if(rstatus==-1) {
      printf("out of grid: %5d %5.3f %6.3f\n",k,x,y);
      }
    for(m=0;m<count;m++)
      G[k+nodes[m]*dimD]=weights[m];
    }

  for(k=0;k<dimD;k++) {
    dum=0.;
    for(m=0;m<dimM;m++)
      dum+=G[m*dimD+k]*mss_cls[m];
    if(fabs(dum-TPdata.data[k].values[24]) > 0.01) {
      printf("%6d %lf %6.3f %6.3f %6.3f %6.3f %3d \n",k,TPdata.data[k].time,TPdata.data[k].lon,TPdata.data[k].lat,TPdata.data[k].values[24],dum,TPdata.data[k].cycle);
      }
    }

/*-----------------------------------------------------------------------------
  compute G transpose, i.e. Gx */

  rstatus=gmatrix_transpose(G, tG, dimD, dimM);


/*-----------------------------------------------------------------------------
  simple test : set diagonal covariance matrix Cm
  goto next;
*/


/*-----------------------------------------------------------------------------
  simple test : set diagonal covariance matrix Cm*/

  rstatus=set_byderivative(Cm,  dimM,  mssgrid);
  goto next;

  rstatus=set_gaussian(Cm,  dimM,  mssgrid);
  goto next;

  rstatus=set_topographic(Cm,  dimM,  mssgrid, mss_cls);
  goto next;

 next:

  printf("model covariance matrix done...\n");

/*-----------------------------------------------------------------------------
  internal consistency chek, Cm must be symmetric */
  rstatus=dgmatrix_symmetry(Cm, dimM);

  printf("rank of Cm     = %d\n",dimM);
  printf("upper bw of Cm = %d\n",dgmatrix_ubw(Cm,dimM,dimM));
  printf("lower bw of Cm = %d\n",dgmatrix_lbw(Cm,dimM,dimM));


  //CETTE BOUCLE NE SERT A RIEN CAR  VARIANCE1 N'EST PAS UTILISEE
  variance1 =(float *)malloc(dimM*sizeof(float));
  for(k=0;k<dimM;k++) variance1[k]=Cm[k*dimM+k];
  free(variance1);variance1=NULL;

/*-----------------------------------------------------------------------------
  compute Gx x Cddeg */

  tmp1 =(float *)malloc(dimM*dimD*sizeof(float));
  for(k=0;k<dimM*dimD;k++) tmp1[k]=0;

  for(l=0;l<dimD;l++) {
    for(k=0;k<dimM;k++) {
      for(m=l-diagD;m<=l+diagD;m++) {
        n=l-m+diagD;                    /* m,l */
/*      tmp1[l*dimM+k]+=tG[m*dimM+k]*Cd[l*dimD+m]; obsolete */
        tmp1[l*dimM+k]+=tG[m*dimM+k]*Cd[l*ldaD+n];
        }
      }
    }
  rstatus=gmatrix_ubw(tG,dimM,dimD);
  rstatus=gmatrix_ubw(tmp1,dimM,dimD);
  free(tG);tG=NULL;

/*-----------------------------------------------------------------------------
  compute Gx x Cddeg x G */

  tmp2 =(float *)malloc(dimM*dimM*sizeof(float));
  for(k=0;k<dimM*dimM;k++) tmp2[k]=0;

  svector =(float *)malloc(dimD*sizeof(float));
  index  =(int *)malloc(dimD*sizeof(int));

  for(k=0;k<dimM;k++) {
/*-----------------------------------------------------------------------------
    get row k of tmp1 */
    nindex=0;
    for(m=0;m<dimD;m++) {
      svector[m]=tmp1[m*dimM+k];
      if(svector[m] !=0. ) {
        nindex++;
        index[nindex-1]=m;
        }
      }
    if(nindex==0) {
      continue;
      }

    if(nindex < 100) {
      for(l=0;l<dimM;l++) {
/*-----------------------------------------------------------------------------
      scalar product with col l of G */
        for (n=0;n<nindex;n++) {
          m=index[n];
          tmp2[l*dimM+k]+= svector[m]*G[l*dimD+m];
          }
        }
      }
    else {
      for(l=0;l<dimM;l++)

#if ATLAS == 1
   tmp2[l*dimM+k]=cblas_dsdot(dimD, svector, 1,&(G[l*dimD]) , 1);
#elif CBLAS == 1
   tmp2[l*dimM+k]=cblas_dsdot(dimD, svector, 1,&(G[l*dimD]) , 1);
#else
#error "stuff missing"
#endif
      }
    }

/*-----------------------------------------------------------------------------
  old method
  for(k=0;k<dimM;k++) {
    for(m=0;m<dimD;m++) vector[m]=tmp1[m*dimM+k];
    norm= sdot( dimD, vector, 1, vector, 1);
    if(norm==0.) {
      continue;
      }
    for(l=0;l<dimM;l++) {
      tmp2[l*dimM+k]= sdot( dimD, vector, 1, &(G[l*dimD]), 1);
      }
    }
*/
  free(svector);svector=NULL;
  free(index);index=NULL;

/*-----------------------------------------------------------------------------
  internal consistency chek, tmp2 must be symmetric */
  rstatus=sgmatrix_symmetry(tmp2, dimM);

  printf("upper bw of  Gx x Cddeg x G %d\n",gmatrix_ubw(tmp2,dimM,dimM));
  printf("lower bw of  Gx x Cddeg x G %d\n",gmatrix_lbw(tmp2,dimM,dimM));

/*-----------------------------------------------------------------------------
  compute Gx x Cddeg x G + Cmdeg */

  A =(double *)malloc(dimM*dimM*sizeof(double));
  for(k=0;k<dimM*dimM;k++) A[k]=tmp2[k]+Cm[k];
/*   for(k=0;k<dimM*dimM;k++) A[k]=tmp2[k]; */
/*   for(k=0;k<dimM*dimM;k++) A[k]=Cm[k]; */
  free(tmp2);tmp2=NULL;

/*-----------------------------------------------------------------------------
  internal consistency chek, tmp2 must be symmetric */
  rstatus=dgmatrix_symmetry(A, dimM);

/*-----------------------------------------------------------------------------
  compute d-Gm  */

  model=mss_cls;

  observation=(float *)malloc(dimD*sizeof(float));
  for(k=0;k<dimD;k++) {
/*     observation[k]=TPdata.data[k].values[0]-TPdata.data[k].values[1]-TPdata.data[k].values[12]; LR, 07/02/05, zone AMSUD!!! */
    observation[k]=TPdata.data[k].values[0]-TPdata.data[k].values[1]-TPdata.data[k].values[18];
    observation[k]=TPdata.data[k].values[29];
    }

  for(k=0;k<dimD;k++) {
      error[k]=observation[k];
    dum=0.;
/*-----------------------------------------------------------------------------
    compute row k, col m  of G x model  */
    for(m=0;m<dimM;m++) dum+=G[m*dimD+k]*model[m];
    error[k]-=dum;
    }

#if ATLAS == 1
 dnorm= ( cblas_ddot(dimD, error, 1,error, 1) / dimD );
#elif CBLAS == 1
 dnorm= ( cblas_ddot(dimD, error, 1,error, 1) / dimD );
#else
#error "stuff missing"
#endif

  printf("prior error euclidian norm (rms) : %f \n",sqrt(dnorm));

/*-----------------------------------------------------------------------------
  compute Gx x Cddeg x (d-Gm)  */

  b =(double *)malloc(dimM*sizeof(double));
  for(k=0;k<dimM;k++) {
    b[k]=0.;
/*-----------------------------------------------------------------------------
    compute row k, col m  of (Gx x Cddeg) x row m of error */
    for(m=0;m<dimD;m++) b[k]+=tmp1[m*dimM+k]*error[m];
    }
  free(tmp1);tmp1=NULL;

/*-----------------------------------------------------------------------------
  compute prior J */

  dvector =(double *)malloc(dimD*sizeof(double));
  for(k=0;k<dimD;k++) {
    dvector[k]=0;
    for(m=k-diagD;m<=k+diagD;m++) {
      n=k-m+diagD;
      dvector[k]+=Cd[k*ldaD+n]*error[m];
      }
    }

#if ATLAS == 1
 dnorm=cblas_ddot(dimD, dvector, 1,error, 1);
#elif CBLAS == 1
 dnorm=cblas_ddot(dimD, dvector, 1,error, 1);
#else
   #error
#endif
     free(dvector);dvector=NULL;

  printf("prior error penalty function : %f %f \n",dnorm,dnorm/dimD);

/*-----------------------------------------------------------------------------
  convert linear problem matrix into band matrix shape*/

/*-----------------------------------------------------------------------------
  solve the linear problem using ??? later should be Cholesky*/

  job=0;

/*   printf("upper bw of A %d\n",gmatrix_ubw(A,dimM,dimM)); */
/*   printf("lower bw of A %d\n",gmatrix_lbw(A,dimM,dimM)); */

  neq=dimM;
  pivot =(int *)malloc(dimM*sizeof(int));

  rstatus=poc_getrf(neq,A,pivot);
  int nrhs=1;
  rstatus=poc_getrs(neq,nrhs,A,pivot,b);

  free(A);A=NULL;
  free(pivot);pivot=NULL;


/*-----------------------------------------------------------------------------
  compute analysis  */

  for(k=0;k<dimM;k++) {
    innovation[k]=b[k];
    analysis[k]=model[k]+b[k];
    }

/*-----------------------------------------------------------------------------
  compute posterior d-Gm  */

  model=analysis;

  for(k=0;k<dimD;k++) {
    error[k]=observation[k];
    dum=0.;
/*-----------------------------------------------------------------------------
    compute row k, col m  of G x model  */
    for(m=0;m<dimM;m++) dum+=G[m*dimD+k]*model[m];
    error[k]-=dum;
    }
  free(observation);observation=NULL;

#if ATLAS == 1
 dnorm=cblas_ddot(dimD, error, 1,error, 1) / dimD;
#elif CBLAS == 1
 dnorm=cblas_ddot(dimD, error, 1,error, 1) / dimD;
#else
#error "stuff missing"
#endif
  printf("posterior error euclidian norm (rms) : %f \n",sqrt(dnorm));

/*-----------------------------------------------------------------------------
  compute posterior J */

/*-----------------------------------------------------------------------------
  data departure = (d-Gm')x Cddeg (d-Gm') */
  dvector =(double *)malloc(dimD*sizeof(double));
  for(k=0;k<dimD;k++) {
    dvector[k]=0;
    for(m=k-diagD;m<=k+diagD;m++) {
      n=k-m+diagD;
      dvector[k]+=Cd[k*ldaD+n]*error[m];
      }
    }

#if ATLAS == 1
 dnorm=cblas_ddot(dimD,dvector , 1,error, 1);
#elif CBLAS == 1
 dnorm=cblas_ddot(dimD,dvector , 1,error, 1);
#else
#error "stuff missing"
#endif
  free(dvector);dvector=NULL;
   free(Cd);Cd=NULL;
  printf("posterior error penalty function : %f %f \n",dnorm,dnorm/dimD);

/*-----------------------------------------------------------------------------
  model departure = (m0-m')x Cmdeg (m0-m') */
  dvector =(double *)malloc(dimM*sizeof(double));
  for(k=0;k<dimM;k++) {
    dvector[k]=0;
    for(m=0;m<dimM;m++) {
      dvector[k]+=Cm[m*dimM+k]*b[m];
      }
    }

#if ATLAS == 1
 dnorm+=cblas_ddot(dimD,dvector , 1,b, 1);
#elif CBLAS == 1
 dnorm+=cblas_ddot(dimD,dvector , 1,b, 1);
#else
#error "stuff missing"
#endif
  free(dvector);dvector=NULL;
  free(Cm);Cm=NULL;

  printf("posterior error penalty function : %f %f \n",dnorm,dnorm/dimD);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int compute_annual(float* h, double *t, float mask, int count)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,n,status;
  double *time,*serie,*residuals,mean;
  spectrum_t s;
  statistic_t sdum;

/*-----------------------------------------------------------------------------
  remove annual signal*/

  status=spectrum_init(&s);
  status=addwave(&s, wSa);
  spectrum_terminate(s);

  serie    =(double *)malloc(count*sizeof(double));
  residuals=(double *)malloc(count*sizeof(double));
  time=(double *)malloc(count*sizeof(double));

  n=0;
  for(k=0;k<count;k++) {
    if(h[k]!=mask) {
      serie[n]=h[k];
      time[n]=t[k];
      n++;
      }
    }

  status=harmonic_analysis_old(serie,time,residuals,&mean,n, s);
  sdum=get_statistics(serie,     (double) 9999.9, n);
  sdum=get_statistics(residuals, (double) 9999.9, n);

  free(serie);serie=NULL;
  free(residuals);residuals=NULL;
  free(time);time=NULL;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

grid_t get_trackgridsioux(double *lon,double *lat, serie_t data,grid_t *nominal_track,char sens[3])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;
  grid_t cgrid;
  double dt,dp,dn;
  double tx,ty,nx,ny;
  double x,y,*d,*alpha,*dtime;
  double sum1,sum2,a,c,kappa;
  grid_t backbone;
  int i,j,k,dx,dy;
  double t0,p0;


  int zone_initialised=0;

/* cartesian, rotated grid */

  cgrid.ymin=   +0.0;
  cgrid.ymax= +200.0 ;
  cgrid.xmin=   -2.0 ;
  cgrid.xmax=   +2.0 ;
  cgrid.dx  =    1.0;
  cgrid.dy  =    2.0;

  dt=lon[1]-lon[0];
  dp=lat[1]-lat[0];
  dn=sqrt(dt*dt+dp*dp);
  dt/=dn;
  dp/=dn;

  tx=dt;
  ty=dp;

  nx=-dp;
  ny= dt;

  grid.dx  =  1.0/100.;
  grid.nx  =   5;
  dx=2;

  grid.dx  =  0.25/100.;
  grid.nx  =   9;
  dx=4;

  grid.dy  =  0.5/100.;
  dy=(int)( 0.1/grid.dy );
  grid.ny  = (int)( 2*dy+dn/grid.dy );

  grid.dx  =  1.5/100.;
  grid.nx  =   3;
  dx=(int)( NINT(grid.dx-1)/2 );
  dx=1;

  grid.dx  =   1.0/50.;
  grid.nx  =   5;
  dx=2;

  grid.dy  =  1.0/50.;
  dy=(int)( 0.1/grid.dy );
  grid.ny  = (int)( 2*dy+dn/grid.dy );

/* test !!!!!!  */
/*   grid.dx  =  3.0/100.; */
/*   grid.nx  =   3; */
/*   dx=1; */
/*   grid.dy  =  0.5/100.; */
/*   dy=0.1/grid.dy; */
/*   grid.ny  = 2*dy+dn/grid.dy; */

  
  grid.modeH  = 2;

  t0=lon[0]+dx*nx*grid.dx-dy*tx*grid.dy;
  p0=lat[0]+dx*ny*grid.dx-dy*ty*grid.dy;

  grid.x=(double *)malloc(grid.nx*grid.ny*sizeof(double));
  grid.y=(double *)malloc(grid.nx*grid.ny*sizeof(double));

  grid.ymin= +1.e+35;
  grid.ymax= -1.e+35;
  grid.xmin= +1.e+35;
  grid.xmax= -1.e+35;

  for (j=0;j<grid.ny;j++)
    for (i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      grid.x[k]=t0-i*nx*grid.dx+j*tx*grid.dy;
      grid.y[k]=p0-i*ny*grid.dx+j*ty*grid.dy;
      grid.ymin=MIN(grid.ymin,grid.y[k]);
      grid.ymax=MAX(grid.ymax,grid.y[k]);
      grid.xmin=MIN(grid.xmin,grid.x[k]);
      grid.xmax=MAX(grid.xmax,grid.x[k]);
      }

  t0=lon[0]-dy*tx*grid.dy;
  p0=lat[0]-dy*ty*grid.dy;

 d=(double *)malloc(data.count*sizeof(double));
  alpha=(double *)malloc(data.count*sizeof(double));
  for (k=0;k<data.count;k++) {
    x=data.data[k].lon;
    y=data.data[k].lat;
    d[k]=(x-t0)*nx+(y-p0)*ny;
    alpha[k]=((x-t0)*tx+(y-p0)*ty)/dn;
    }

  sum1=0;
  sum2=0;

  for (k=0;k<data.count;k++) {
    c=(1.-alpha[k])*alpha[k];
    sum2+=d[k]*c;
    sum1+=c*c;
    }

  free(alpha);alpha=NULL;
  free(d);d=NULL;

  kappa=sum2/sum1;

  backbone.ny=grid.ny*2;
  backbone.x=(double *)malloc(backbone.ny*sizeof(double));
  backbone.y=(double *)malloc(backbone.ny*sizeof(double));
  dtime=(double *)calloc(backbone.ny,sizeof(double));
  for(i=0;i<backbone.ny;i++)dtime[i]=i/2.0;
 
//    interpTrace( sens, lon[0]*DEG2RAD, lat[0]*DEG2RAD, backbone.ny, dtime, backbone.x, backbone.y );
//    for(i=0;i<backbone.ny;i++){backbone.x[i]*=RAD2DEG;backbone.y[i]*=RAD2DEG;}

   for (j=0;j<backbone.ny;j++) {
      a=(double) j/(double)backbone.ny;
      c=(1.-a)*a;
      backbone.x[j]=t0+j*tx*grid.dy+c*kappa*nx;
      backbone.y[j]=p0+j*ty*grid.dy+c*kappa*ny;
      }

  grid.ymin= +1.e+35;
  grid.ymax= -1.e+35;
  grid.xmin= +1.e+35;
  grid.xmax= -1.e+35;

  for (j=0;j<grid.ny;j++)
    for (i=0;i<grid.nx;i++) {
      k=j*grid.nx+i;
      grid.x[k]=backbone.x[j]-(i-dx)*nx*grid.dx;
      grid.y[k]=backbone.y[j]-(i-dx)*ny*grid.dx;
      grid.ymin=MIN(grid.ymin,grid.y[k]);
      grid.ymax=MAX(grid.ymax,grid.y[k]);
      grid.xmin=MIN(grid.xmin,grid.x[k]);
      grid.xmax=MAX(grid.xmax,grid.x[k]);
      }

  zone_initialised=1;
   
  if(zone_initialised!=1) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

   j= compute_nominal_track(&backbone,nominal_track);

  free(backbone.x);backbone.x=NULL;
  free(backbone.y);backbone.y=NULL;
  free(dtime);dtime=NULL;

  return(grid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void get_extremities(char *sat,int k,serie_t TPdata,double *lon,double *lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count;
  double x,y;

  if (strcasecmp(sat,"GFO")!=0)
    /*-----------------------------------------------------------------------------
      all satellites except GFO */
    {
      if(k%2==0) /*descending*/
        {
          lon[0]=+1.e+35;
          lat[0]=-1.e+35;
          lon[1]=-1.e+35;
          lat[1]=+1.e+35;
          for(k=0;k<TPdata.count;k++) {
              lon[0]=MIN(TPdata.data[k].lon,lon[0]);
              lon[1]=MAX(TPdata.data[k].lon,lon[1]);
              lat[0]=MAX(TPdata.data[k].lat,lat[0]);
              lat[1]=MIN(TPdata.data[k].lat,lat[1]);
            }
        }
      else
        {
          lon[0]=+1.e+35;
          lat[0]=+1.e+35;
          lon[1]=-1.e+35;
          lat[1]=-1.e+35;
          for(k=0;k<TPdata.count;k++) {
              lon[0]=MIN(TPdata.data[k].lon,lon[0]);
              lon[1]=MAX(TPdata.data[k].lon,lon[1]);
              lat[0]=MIN(TPdata.data[k].lat,lat[0]);
              lat[1]=MAX(TPdata.data[k].lat,lat[1]);
            }
        }

    }
  else
    /*-----------------------------------------------------------------------------
      GFO */
    {
      if(k%2==0) /*descending*/
        {
          lon[0]=+1.e+35;
          lat[0]=+1.e+35;
          lon[1]=-1.e+35;
          lat[1]=-1.e+35;
          for(k=0;k<TPdata.count;k++) {
              lon[0]=MIN(TPdata.data[k].lon,lon[0]);
              lon[1]=MAX(TPdata.data[k].lon,lon[1]);
              lat[0]=MIN(TPdata.data[k].lat,lat[0]);
              lat[1]=MAX(TPdata.data[k].lat,lat[1]);
            }
        }
      else
        {
          lon[0]=+1.e+35;
          lat[0]=-1.e+35;
          lon[1]=-1.e+35;
          lat[1]=+1.e+35;
          for(k=0;k<TPdata.count;k++) {
              lon[0]=MIN(TPdata.data[k].lon,lon[0]);
              lon[1]=MAX(TPdata.data[k].lon,lon[1]);
              lat[0]=MAX(TPdata.data[k].lat,lat[0]);
              lat[1]=MIN(TPdata.data[k].lat,lat[1]);
            }
        }

    }

/*   mssgrid=get_trackgrid(lon,lat); */
  printf("extrimities :  %f %f %f %f\n",lon[0],lat[0],lon[1],lat[1]);


  /*------------------------------------------------------------------------------
    smooth extremities computation */


    count=0;
    x=0;
    y=0;
    for(k=0;k<TPdata.count;k++) {
      if((fabs(lon[0]-TPdata.data[k].lon) < 0.05) || (fabs(TPdata.data[k].lat-lat[0]) < 0.05)) {
        x+=TPdata.data[k].lon;
        y+=TPdata.data[k].lat;
        count++;
        }
      }
    lon[0]=x/count;
    lat[0]=y/count;

    count=0;
    x=0;
    y=0;
    for(k=0;k<TPdata.count;k++) {
      if((fabs(lon[1]-TPdata.data[k].lon) < 0.05) || (fabs(TPdata.data[k].lat-lat[1]) < 0.05)) {
        x+=TPdata.data[k].lon;
        y+=TPdata.data[k].lat;
        count++;
        }
      }
    lon[1]=x/count;
    lat[1]=y/count;

/*
  projection=geo_mercator_init(lon[1],lat[1], (double) 0.0, (double) 90.0);
  status=geo_mercator_directe(projection,lon[0],lat[0],&x,&y);
  status=geo_mercator_directe(projection,lon[1],lat[1],&x,&y);
  x=geo_distance(lon[0],lat[0],lon[1],lat[1]);
*/

  printf("extrimities :  %f %f %f %f\n",lon[0],lat[0],lon[1],lat[1]);
}
