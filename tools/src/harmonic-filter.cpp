
/*******************************************************************************

  T-UGO tools, 2006-2019

  Unstructured Ocean Grid initiative

*******************************************************************************/

#include <vector>

#include "tools-structures.h"

#include "maths.h"
#include "tides.h"
#include "lanczos.h"
#include "loess.h"

#include "mgr-converter.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_load_tseries(const char *filename, vector<int> & columns, tseries_t & series, char unit,const char *header_format, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, n,ndata,ncol;/* <index, number of lines, number of time steps */
  char line[300];
  FILE *in=NULL;
  char key[32];
  char * pos;
  double *z;
  
  in=fopen(filename,"r");
  if (in == NULL) TRAP_ERR_EXIT(-1,"unable to open file: %s \n",filename);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  get max number of data/column in file
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  n=0;
  ndata=0;
  while(!feof(in)) {
    fgets(line,sizeof(line),in);
    n++;
    if(line[0]=='#') continue;
    ndata++;
    
    i=count_token(line, " \t");
    updatemax(&ncol,i-1);
    }
  n--;
  ndata--;
  
  rewind(in);
  
  z=new double[ncol];

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  allocate memory
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  series=tseries_t(ndata, columns.size());
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  read data
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  ndata=0;
  for (i=0; i<n;i++) {
    fgets(line,sizeof(line),in);
/*------------------------------------------------------------------------------
    first line may be header */
    if (i==0 && (mooring==NULL || mooring->name==NULL || mooring->code == -1)) {
      if (header_format)
        mgr_parse_mooring_header(header_format, line, mooring);
      else if(line[0]!='#') {
        if (mooring->name==NULL)
          mooring->name=new char[256];
        sscanf(line, "%s %s %d %lf %lf %lf", key, mooring->name, &(mooring->code), &(mooring->lon),&(mooring->lat),&(mooring->depth));
        }
      }
/*------------------------------------------------------------------------------
    standard line */
    if(line[0]!='#') {
/*------------------------------------------------------------------------------
      read time */
      if(line[0]=='\n') continue;
      for(int k=0;k<300;k++) {
        if(line[k]=='\t') line[k]=' ';
        if(line[k]=='\n') break;
        }
      pos = strtok (line," ");
      int nitems=sscanf(pos, " %lf", &series.t[ndata]);
/*------------------------------------------------------------------------------
      loop over data */
      int j = 0;
      do {
        pos = strtok (NULL," ");
        if (pos==0) pos = strtok (NULL,"\t");
        if (pos==0) TRAP_ERR_EXIT(-1,"** ERROR: bad format file: %s \n",filename);
        nitems=sscanf(pos, " %lf", &z[j]);
        j++;
        } while (j<ncol);
      for(int k=0;k<columns.size();k++) series.x[k][ndata]=z[columns[k]];
      ndata++;
      }
    }

  fclose(in);
  
  delete[] z;
  
  return(ndata);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HarmonicFilter_init(double *t, const vector<double> & f, double* &H, double* &LF, size_t ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Tunable, Demerliac type f filtering:
  
  y=Pa P is prediction matrix [cos sin cos sin ...]
  
  min[sum(h-y)簡] -> tPPa=tPh y=[P x inv(tPP) x tP] h
  
  Low-pass filter : LF = [I - P x inv(tPP) x tP ] h
  
  Nota:
  
    if P invertible, lf=0. To eliminate T=2, 3, ..., 8 hours
 
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int  status;
  int verbose=0;
  double *P, *HF=0;
//   double *H=0;
  vector<double> omega;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  initialise prediction matrix

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  for(int k=0; k<f.size(); k++) {
    double frequency=f[k];
    omega.push_back(frequency*2.0*M_PI);
    }
  
  size_t ncols=2*f.size();
  size_t nrows=ndata;
  
  P=new double[ncols*nrows];
  
  for(int n=0;n<ndata;n++) {
    double time=t[n];
    for(int k=0; k<f.size(); k++) {
      int col;
      double phase=fmod(omega[k]*time, 2.0*M_PI);
      col=2*k;
//       P[n*ncols+col]=cos(phase);
      P[col*nrows+n]=cos(phase);  /* fortran ordering  : row n, column col */
      col=2*k+1;
//       P[n*ncols+col]=sin(phase);
      P[col*nrows+n]=sin(phase);   /* fortran ordering  : row n, column col */
     
      }
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  initialise harmonic matrix

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status=matrix_least_square(P, nrows, ncols, H, verbose);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  initialise high-pass filter matrix

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  status=matrix_product(P, H, nrows, ncols, ncols, nrows, HF, verbose);

  LF=new double[nrows*nrows];
  
  for(int n=0;n<nrows;n++) {
    double w=0;
    for(int m=0;m<nrows;m++) {
      LF[m*nrows+n]=-HF[m*nrows+n];
      w+=HF[m*nrows+n];
      }
//     printf("%d %lf\n",n,w);
/*------------------------------------------------------------------------------
    diagonal term */
    LF[n*nrows+n]+=1.0;
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HarmonicFilter(double *t, double *buffer, double *filtered, double mask, size_t ndata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  Tunable, Demerliac type filtering:
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int  status;
  double *LF=0, *H=0;
  vector<double> f;
  int width, hw;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  check for harmonic-specified input

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define CHECKS 0

#if CHECKS
  int verbose=0;
  
  for(int n=0;n<ndata;n++) {
    double T;
    T=12. / 24.;
    buffer[n]=0.5*cos(2.*M_PI*t[n]/T);
    T=12.47 / 24.;
    buffer[n]+=1.0*cos(2.*M_PI*t[n]/T+1.0);
//     T=24. / 24.;
//     buffer[n]+=1.0*cos(2.*M_PI*t[n]/T);
    }
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  select frequencies (cycles/julian days)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
//   for(int k=4; k<13; k++) {
//     double T=(double) k / 24.;
//     double frequency=1.0/T;
//     f.push_back(frequency);
//     }
  
  double T, frequency;
  
//   T=2.18181818181818181818/24.; frequency=1.0/T;
//   f.push_back(frequency);
  
//   T=2.4/24.; frequency=1.0/T;
//   f.push_back(frequency);
  
//   T=2.666666666/24.; frequency=1.0/T;
//   f.push_back(frequency);
  
#if HarmonicFilterSpectrumType == 1
  T=3.42857/24.; frequency=1.0/T;
  f.push_back(frequency);
  
  T=4.8/24.; frequency=1.0/T;
  f.push_back(frequency);
  
  T=6./24.; frequency=1.0/T;
  f.push_back(frequency);
  
  T=8./24.; frequency=1.0/T;
  f.push_back(frequency);
  
  T=12./24.; frequency=1.0/T;
  f.push_back(frequency);
  
  T=12.5/24.; frequency=1.0/T;
  f.push_back(frequency);
  
  T=13/24.; frequency=1.0/T;
  f.push_back(frequency);
  
  T=1.0; frequency=1.0/T;
  f.push_back(frequency);
#else
  spectrum_t WaveList;
  
#if HarmonicFilterSpectrumType == 2
  char *waves[]={"S1","M2","S2","S3","S4","S5","S6","S7","S8","S9","S10","S11",0};
#elif HarmonicFilterSpectrumType == 3
  char *waves[]={"S1","M2","S2","S3","S4","S5","S7","S9",0};
#elif HarmonicFilterSpectrumType == 4
  char *waves[]={"S1","M2","S2","S3","S4","S5",0};
#endif
  
  WaveList.init(initialize_tide(),waves);
  for(int k=0;k<WaveList.n;k++){
    frequency=WaveList.waves[k].omega/15.;
    f.push_back(frequency);
    }
#endif
  
  for(int k=0;k<f.size();k++) printf("%s : %d %lf cycles/day\n", __func__, k, f[k]);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  init filter
  
  filter window will also dictate separation criterion (frequency precision)
  
  

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  filtering window (in days) */
  double window=3.0;
  
  double sampling=get_timesampling(t, mask, ndata);
  
  hw=nint(window/sampling);
  width=2*hw+1;
  
  status=HarmonicFilter_init(t, f, H, LF, width);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  check for harmonic-specified input

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if CHECKS
  a=new double[2*f.size()];
  
  status=matrix_operation(H, buffer, (int) 2*f.size(), width, a, verbose);
  
  for(int k=0;k<f.size();k++) printf("%s : %d %lf\n", __func__, k, sqrt(a[2*k]*a[2*k]+a[2*k+1]*a[2*k+1]));
  
  status=matrix_operation(H, buffer+1024, (int) 2*f.size(), width, a, verbose);
  
  for(int k=0;k<f.size();k++) printf("%s : %d %lf\n", __func__, k, sqrt(a[2*k]*a[2*k]+a[2*k+1]*a[2*k+1]));
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  process filtering

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  for(int n=0;n<ndata;n++) {
    double w=0;
    filtered[n]=0.0;
    for(int k=0;k<width;k++) {
      int m=n+k-hw;
      if(m<0)      continue;
      if(m>=ndata) continue;
      filtered[n]+=buffer[m]*LF[hw*width+k];
      w+=LF[hw*width+k];
      }
//     printf("%d %lf\n",n,w);
//     filtered[n]/=w;
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HarmonicFilter_experiment()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  ERA5 experiment:
  
  from december 2018 checks, it seems that atmospheric pressure high frequency
  is aliased in the hourly archives
  
  experimental fix:
  
  -keep Sk harmonics as they seem to be meaningfull (via Demerliac)
  
  -remove continuum VHF signal (T<6h, via Lanczos)
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int  status;
  vector<double> f;
  
  vector<int> columns;
  tseries_t series;
  mooring_t mooring;
  string format, time_format;
  string rootname;
  bool gnu_nice=true;
  
  status=mgr_init_formats();
  
  columns.push_back(0);
  status=mgr_load_tseries("sample.15", columns, series, 'm', 0, &mooring);

  series.mask=1.e+10;

/*------------------------------------------------------------------------------
  add two parameter */  
  series.redim(4);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  first isolate HF signal with 8h period cut-off
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  double sampling=get_timesampling(series.t, series.mask, series.n);
  double *weightBF,half;
  int nWeightBF;
 
  half=8.0/24./2.0;
  lanczos1D_init(sampling,series.n,half,&weightBF,&nWeightBF);
  
  Loess1D_BF_nomask(series.x[0],series.n,weightBF,nWeightBF,series.x[1]);
  /* 1 : Lanczos LF */
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  then isolate harmonic signal (S1, S2, S3, ...) in HF
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  for(int n=0;n<series.n;n++) {
    series.x[2][n]=series.x[0][n]-series.x[1][n];
    }
  /* 2 : Lanczos HF */
  
  status=HarmonicFilter(series.t, series.x[2], series.x[3], series.mask, series.n);
  /* 3 : continuum of Lanczos HF */
  
  for(int n=0;n<series.n;n++) {
    series.x[4][n]=series.x[3][n];
    series.x[3][n]=series.x[2][n]-series.x[3][n];
    }

  range_t<size_t> range(6*nWeightBF, series.n-6*nWeightBF);
  statistic_t s=get_statistics(series.x[4],series.mask,range,1);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  finally put harmonic signal back in BF
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  for(int n=0;n<series.n;n++) {
    series.x[1][n]=series.x[1][n]+series.x[3][n];
    }
  
  /*
  0 : original
  1 : Lanczos LF + harmonic of Lanczos HF
  2 : Lanczos HF
  3 : harmonic of Lanczos HF
  4 : continuum of Lanczos HF
  */
  
  format=FORMAT_NAME_GNU;
  rootname="filtered";
  status=mgr_save_timeserie(rootname, mooring, series, 'm', format, gnu_nice, time_format);
  
  series.destroy();
  
  return(status);
}
