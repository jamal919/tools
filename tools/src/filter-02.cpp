/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"

#include "tools-define.h"
#include "poc-time.h"
#include "filter.h"
#include "tides.h"
#include "statistic.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_coef(const char *filter, double **coef, int *ncoef)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  int i;
  double s;

  if(strcmp(filter,"demerliac") == 0 ) {
    *coef= new double[71];
/* a(-k) = a(k) for the Demerliac coef */
/* here, the coef are stored from coef[0] -> a(-35) to coef[70] -> a(35) */
    (*coef)[0]=1;
    (*coef)[1]=3;
    (*coef)[2]=8;
    (*coef)[3]=15;
    (*coef)[4]=21;
    (*coef)[5]=32;
    (*coef)[6]=45;
    (*coef)[7]=55;
    (*coef)[8]=72;
    (*coef)[9]=91;
    (*coef)[10]=105;
    (*coef)[11]=128;
    (*coef)[12]=153;
    (*coef)[13]=171;
    (*coef)[14]=200;
    (*coef)[15]=231;
    (*coef)[16]=253;
    (*coef)[17]=288;
    (*coef)[18]=325;
    (*coef)[19]=351;
    (*coef)[20]=392;
    (*coef)[21]=435;
    (*coef)[22]=465;
    (*coef)[23]=512;
    (*coef)[24]=558;
    (*coef)[25]=586;
    (*coef)[26]=624;
    (*coef)[27]=658;
    (*coef)[28]=678;
    (*coef)[29]=704;
    (*coef)[30]=726;
    (*coef)[31]=738;
    (*coef)[32]=752;
    (*coef)[33]=762;
    (*coef)[34]=766;
    (*coef)[35]=768;

    for(i=1;i<36;i++) (*coef)[35+i]=(*coef)[35-i];

    s=0;
    for (i=0;i<71;i++) s+=(*coef)[i];

    s=24576;
    for (i=0;i<71;i++) (*coef)[i]/=s;
    *ncoef=35;
    status=0;
    }

  if(strcmp(filter,"godin") == 0 ) {
    *coef= new double[71];
/* a(-k) = a(k) for the Godin coef */
/* here, the coef are stored from coef[0] -> a(-35) to coef[70] -> a(35) */
    (*coef)[0]=1;
    (*coef)[1]=3;
    (*coef)[2]=6;
    (*coef)[3]=10;
    (*coef)[4]=15;
    (*coef)[5]=21;
    (*coef)[6]=28;
    (*coef)[7]=36;
    (*coef)[8]=45;
    (*coef)[9]=55;
    (*coef)[10]=66;
    (*coef)[11]=78;
    (*coef)[12]=91;
    (*coef)[13]=105;
    (*coef)[14]=120;
    (*coef)[15]=136;
    (*coef)[16]=153;
    (*coef)[17]=171;
    (*coef)[18]=190;
    (*coef)[19]=210;
    (*coef)[20]=231;
    (*coef)[21]=253;
    (*coef)[22]=276;
    (*coef)[23]=300;
    (*coef)[24]=323;
    (*coef)[25]=344;
    (*coef)[26]=363;
    (*coef)[27]=380;
    (*coef)[28]=395;
    (*coef)[29]=408;
    (*coef)[30]=419;
    (*coef)[31]=428;
    (*coef)[32]=435;
    (*coef)[33]=440;
    (*coef)[34]=443;
    (*coef)[35]=444;

    for(i=1;i<36;i++) (*coef)[35+i]=(*coef)[35-i];

    s=0;
    for (i=0;i<71;i++) s+=(*coef)[i];

    s=14400;
    for (i=0;i<71;i++) (*coef)[i]/=s;

    *ncoef=35;
    status=0;
    }

  if(strcmp(filter,"munk") == 0 ) {
    *coef= new double[49];

/* a(-k) = a(k) for the Munk coef */
/* here, the coef are stored from coef[0] -> a(-24) to coef[48] -> a(24) */
    (*coef)[0] = 13307;
    (*coef)[1] = 30073;
    (*coef)[2] = 47028;
    (*coef)[3] = 60772;
    (*coef)[4] = 72261;
    (*coef)[5] = 85349;
    (*coef)[6] =101603;
    (*coef)[7] =122665;
    (*coef)[8] =146225;
    (*coef)[9] =165525;
    (*coef)[10]=180727;
    (*coef)[11]=195528;
    (*coef)[12]=208050;
    (*coef)[13]=219260;
    (*coef)[14]=234033;
    (*coef)[15]=251492;
    (*coef)[16]=278167;
    (*coef)[17]=300054;
    (*coef)[18]=314959;
    (*coef)[19]=325633;
    (*coef)[20]=338603;
    (*coef)[21]=354118;
    (*coef)[22]=370094;
    (*coef)[23]=386839;
    (*coef)[24]=395287;

    for(i=1;i<25;i++) (*coef)[24+i]=(*coef)[24-i];

    s=0;
    for (i=0;i<49;i++) s+=(*coef)[i];

    s=10000000;/*10e7*/
    for (i=0;i<49;i++) (*coef)[i]/=s;

    *ncoef=24;
    status=0;
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int init_TimeFilter(const char *filtername, tseries_t & filter, double dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int n,half,ncoef;
  double *coefficients, sum;
  tseries_t hourly;
  
  status=init_coef(filtername, &coefficients, &half);
  
  hourly=tseries_t(2*half+1,1);
  
  ncoef=2*half+1;
  
  for(n=0;n<ncoef;n++) {
    hourly.x[0][n]=coefficients[n];
    hourly.t[n]=(double) (n-half) / 24.0;
    }
    
  filter=hourly.resample(dt,1.e+10);
  hourly.destroy();
  
  sum=0;
  for(n=0;n<filter.n;n++) {
    sum+=filter.x[0][n];
    }

  for(n=0;n<filter.n;n++) {
    filter.x[0][n]/=sum;
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *tide_filters(const char *filtername, double *h, double mask, int nvalues )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  designed for hurly timeseries, to be improved */
  int i,j,status;
  int ncoef = 0;
  double *coef_filter = 0;

  double *filtered=new double[nvalues];

  status=init_coef(filtername, &coef_filter, &ncoef);/* !!! shift between the coef and the storage index!!!*/

  i=0;
  while(h[i] == mask) {
    i++;
    } ;
  int start = i;

  for(i = 0; i < nvalues; i++){
    filtered[i] = mask;
    }

  j = 0;
  for(i = start + ncoef; i < nvalues - ncoef; i++) {
    filtered[i] = 0;
    for(j = -ncoef; j < ncoef; j++){
      if(h[i-j] != mask){
        filtered[i] += (coef_filter[j + 35] * h[i - j]);
        }
      else{
        filtered[i]=mask;
        break;
        }
      }
    }

  delete[] coef_filter;
  return(filtered);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *tide_filters(const char *filtername, double dt, double *h, double mask, int nvalues )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  designed for hourly timeseries, to be improved */
  int i,j,half,status;
  double *coef_filter = 0;
  
  tseries_t filter;
  filter.nparam=1;
  
  double *filtered=new double[nvalues];
  
  status=init_TimeFilter(filtername, filter, dt);

  i=0;
  while(h[i] == mask) {
    i++;
    } ;
  int start = i;

  for(i = 0; i < nvalues; i++){
    filtered[i] = mask;
    }

  half=1+filter.n/2;
  
  j = 0;
  for(i = start + half; i < nvalues - half; i++) {
    filtered[i] = 0;
    for(j = -half; j < half+1; j++){
      if(h[i-j] != mask){
        filtered[i] += (filter.x[0][j + half] * h[i - j]);
        }
      else{
        filtered[i]=mask;
        break;
        }
      }
    }

  filter.destroy();
    
  delete[] coef_filter;
  return(filtered);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fill_gaps(double *t, double *h, double mask, int n, double w, bool fit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nLF, status;
  double window;
  double *buffer=new double[n], *LF=new double[n];
  statistic_t   s;
  
/*------------------------------------------------------------------------------
  smooth residuals */
  window=0.125;
  status=loess1d_irregular(n, window, h, t, mask, buffer);
  for (int i=0; i< n; i++) {
    if(fit and h[i] != mask) buffer[i]=h[i];
    }
    
  nLF=0;
  for (int i=0; i< n; i++) {
    if(buffer[i] == mask) nLF++;
    }
    
  window=0.;
  
  while(nLF!=0) {
    window+=w;
    status=loess1d_irregular(n, window, buffer, t, mask, LF);
    for (int i=0; i< n; i++) {
      if(buffer[i] == mask) buffer[i]=LF[i];
      }
    nLF=0;
    for (int i=0; i< n; i++) {
      if(buffer[i] == mask) nLF++;
      }
    }
    
  for (int i=0; i< n; i++) {
    h[i]=buffer[i];
    }
    
  delete[] LF;
  delete[] buffer;

  
  return(window);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int remove_outlayers(double *t, double *h, double mask, double sampling, int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nLF, noutlayers, status;
  double window, delta;
  double *buffer=new double[n], *LF=new double[n], *HF=new double[n];
  double *detided;
  statistic_t   s;
  int previous=0,previous2=0,total;
  int loop=0;
  bool debug=false;
  int verbose=(debug==true);

  total=0;
  loop=0;
      
  redo_peaks:
  
/*------------------------------------------------------------------------------
  fill gaps with low-frequency signal */
  for (int i=0; i< n; i++) {
    buffer[i]=h[i];
    }
  window=fill_gaps(t, buffer, mask, n, 1.0, true);
//   window=0.125;
//   status=loess1d_irregular(n, window, h, t, mask, buffer);
//   for (int i=0; i< n; i++) {
//     if(h[i] != mask) buffer[i]=h[i];
//     }
//     
//   nLF=0;
//   for (int i=0; i< n; i++) {
//     if(buffer[i] == mask) nLF++;
//     }
//     
// /*------------------------------------------------------------------------------
//   */
//   window=0.;
//   while(nLF!=0) {
//     window+=1;
//     status=loess1d_irregular(n, window, buffer, t, mask, LF);
//     for (int i=0; i< n; i++) {
//       if(buffer[i] == mask) buffer[i]=LF[i];
//       }
//     nLF=0;
//     for (int i=0; i< n; i++) {
//       if(buffer[i] == mask) nLF++;
//       }
//     }

/*------------------------------------------------------------------------------
  try Demerliac filter to isolate residual tidal signal */
  detided = tide_filters("demerliac", sampling, buffer, mask, n);
  for (int i=0; i< n; i++) {
    if(detided[i] != mask) {
      HF[i]=buffer[i]-detided[i];
      }
    else {
      HF[i]=mask;
      }
    }
    
  s=get_statistics(HF,mask,n,verbose);
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  remove single peaks 

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  noutlayers=0;
 
  for (int i=1; i< n-1; i++) {
    if((HF[i-1] != mask) && (HF[i] != mask) && (HF[i+1] != mask) ){
      double a=HF[i]-HF[i-1];
      double b=HF[i+1]-HF[i];
      if(a*b<0) {
        if( (fabs(a)>0.15) && (fabs(b)>0.15) ) {
          h[i]=mask;
          noutlayers++;
          printf("%s, peaks : %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf (std=%lf)\n", poctime_sdate_cnes(t[i],'d'), HF[i-2], HF[i-1], HF[i], HF[i+1], HF[i+2], s.std);
          }
        }
      }
    }
  
  if(debug) printf("outlayers (peaks) detection loop %d : #outlayers found %d (previously %d)\n", loop, noutlayers, previous);
  if(noutlayers!=previous) {
    previous=noutlayers;
    total+=noutlayers;
    loop++;
    goto redo_peaks;
    }
  printf("outlayers (step) detection: #loops=%d #outlayers found %d\n", loop, total);
    
  total=0;
  loop=0;
      
  redo_steps:
  
/*------------------------------------------------------------------------------
  fill gaps with low-frequency signal */
  for (int i=0; i< n; i++) {
    buffer[i]=h[i];
    }
  window=fill_gaps(t, buffer, mask, n, 1.0, true);
  
/*------------------------------------------------------------------------------
  try Demerliac filter to isolate residual tidal signal */
  detided = tide_filters("demerliac", sampling, buffer, mask, n);
  for (int i=0; i< n; i++) {
    if(detided[i] != mask and buffer[i] != mask) {
      HF[i]=buffer[i]-detided[i];
      }
    else {
      HF[i]=mask;
      }
    }
    
  s=get_statistics(HF,mask,n,verbose);
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  remove steps
  
  warning: false outlayer detection in estuarine timeseries (falling/rising 
           instant)
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  noutlayers=0;
 
  for (int i=1; i< n; i++) {
    if((HF[i-1] != mask) && (HF[i] != mask) ) {
      double a=HF[i]-HF[i-1];
      double b=HF[i+1]-HF[i];
      if( fabs(a)>3*s.std and fabs(b)>3*s.std) {
//       if( (fabs(a)>0.2)) {
        h[i-1]=mask;
        h[i]  =mask;
//         date_t date=cnesdate(t[i],'d');
        printf("%s, steps : %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf (std=%lf)\n", poctime_sdate_cnes(t[i],'d'), HF[i-2], HF[i-1], HF[i], HF[i+1], HF[i+2], s.std);
        noutlayers++;
        noutlayers++;
        }
      }
    }
  
  if(debug) printf("outlayers (step) detection loop %d : #outlayers found %d (previously %d)\n", loop, noutlayers, previous2);
  if(noutlayers!=previous2) {
    previous2=noutlayers;
    total+=noutlayers;
    loop++;
    goto redo_steps;
    }
    
//   return(noutlayers);
  printf("outlayers (step) detection: #loops=%d #outlayers found %d\n", loop, total);
    
/*------------------------------------------------------------------------------
  remove outliers */
  window=3.;
  status=loess1d_irregular(n, window, HF, t, mask, buffer);
  
  for (int i=0; i< n; i++) {
    if(HF[i] != mask) {
      HF[i]=HF[i]-buffer[i];
      }
    }
  s=get_statistics(HF,mask,n);
  
//   status=mgr_save_timeserie("HF.gnu", 0, *reducedkj, 'm', format_str, gnu_nice);  
  
  total=0;
  noutlayers=0;
  for (int i=0; i< n; i++) {
    if(HF[i] != mask) {
      delta=fabs(HF[i]);
      if(delta>4.*s.std) {
        h[i]=mask;
        noutlayers++;
        printf("%s, outlayers : %6.3lf %6.3lf %6.3lf %6.3lf %6.3lf (std=%lf)\n", poctime_sdate_cnes(t[i],'d'), HF[i-2], HF[i-1], HF[i], HF[i+1], HF[i+2], s.std);
        }
      }
    }
    
  printf("outlayers (step) detection: #loops=%d #outlayers found %d\n", loop, noutlayers);
    
  delete[] HF;
  delete[] LF;
  delete[] buffer;

  return(noutlayers);
 
}
