
/*******************************************************************************

  T-UGO tools, 2006-2016

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Filters/detides TUGOm sample files.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main <a href=#func-members>function</a> for how this works
and the print_help <a href=#func-members>function</a> for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "filter.h"
#include "mgr.h"
#include "functions.h"
#include "tides.def"
#include "spectrum.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : as[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Filters/detides TUGOm sample files.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help  Show this help and exit.\n"
    "  -start  followed by the start date. Default: 01/1900. See DATE FORMATS below\n"
    "  -end  followed by the end date. Default: 12/2100. See DATE FORMATS below\n"
    "  --no-fft  disable costly production of FFTs\n"
    "  -m  followed by samples directory. Default: ./\n"
    "  -i  followed by stations list\n"
    "  -f  followed by timeseries output format: classic (default) or matroos\n"
    "  -o  followed by timeseries output directory. Default: ./\n"
    "  -c  without any effect\n"
    "  -t  produce filtered time series\n"
    "  -d  enable detiding\n"
    "  -r  followed by a wave for the custom spectrum. Default spectrum: S2 .\n"
    "      If the wave can not be fund in an uninitialised list, i.e. always., crash.\n"
    "  -s  followed by a spectrum file\n"
    "\n"
    "BUGS"
    "  See options -c and -r above.\n"
    );
  print_poctime_scan_date_help(0);
  /** \endcode */
}


/*------------------------------------------------------------------------------

  It reads the sample file from the MOG2D model and filter/detide the time
  series

  input:  sample.# where # is the index of the sampling point
          extract.input where the list position/ name is kept
          (presumedly in /export/bartok4/medsea/extract.input)

  output: XXX.h and XXX.p where XXX is the sampling point name

-----------------------------------------------------------------------*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dT[4],dTsec,t0sec,separation;

  double  zi,zr,mean[10],rms[10],count;
  double  **a,**G,a1,p1,a2,p2,d;
  int i,j,k,l,n,fmt=0,status,tg,nbwave;
  int /*nr[4],*/nst,nset=0,test=0,detiding=0,compare=0;
  FILE *out;
  char name[256],tmp[256], datafile[256];
  const char *model=0, *input, *path=NULL,*mgr=NULL,*sfmt=0,*wave,*spectrumFile=0;
  pressure_station *sample;

  spectrum_t spectrum,reference;
  pressure_station *setptr=0;

  harmonic_t x;
  mooring_t mooring;
  int sysinit=0;
  date_t start_date,end_date;
//   double start,end;
  tseries_t serie(4);
  bool doFourier=true;

  fct_echo( argc, argv);

  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);

  spectrum.waves=new tidal_wave[100];
  nbwave=0;

  start_date.init(1900,1);
  end_date.init(2100,12);
  
  n=1;
  while (n < argc) {
    const char *keyword=argv[n];
    
    if(strcmp(keyword,"--help")==0 || strcmp(keyword,"-h")==0) {
      print_help(argv[0]);
      exit(0);
      }
    
    if(strcmp(keyword,"-start") ==0) {
      status=poctime_scan_date(argv[n+1],&start_date,0);
      n++;
      n++;
      continue;
      }
    
    if(strcmp(keyword,"-end") == 0) {
      status=poctime_scan_date(argv[n+1],0,&end_date);
      n++;
      n++;
      continue;
      }
    
    if(strcmp(keyword,"--no-fft") == 0) {
      doFourier=false;
      n++;
      continue;
      }
    
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'm' :
          model= argv[n+1];
          n++;
          n++;
          break;

        case 'i' :
          input= argv[n+1];
          n++;
          n++;
          break;

        case 'f' :
          sfmt= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          path= argv[n+1];
          n++;
          n++;
          break;

        case 'c' :
          mgr= argv[n+1];
          n++;
          n++;
          break;

        case 't' :
          test=1;
          n++;
          break;

        case 'd' :
          detiding=1;
          n++;
          break;

        case 'r' :
          int index;
          detiding=1;
          wave= argv[n+1];
          index = reference.wave_index(wave);
          if (index<0) {
            STDOUT_BASE_LINE("unknown wave index for %s\n", wave);
            exit(-1);
            }
          spectrum.waves[nbwave++] = reference.waves[index];
          n++;
          n++;
          break;

        case 's' :
          spectrumFile= argv[n+1];
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
          break;
        }
        break;

      default:
        STDOUT_BASE_LINE("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
    
    }

  if(path==NULL) {
    path=".";
    printf("use <.> as pathname for default output directorys\n");
    }

  if(mgr!=0) compare=1;

  if(detiding==0) {
    compare=0;
    path=".";
    printf("no deting performed\n");
    }

  if(sfmt==NULL) {
    sfmt="classic";
    printf("use classic for default output format\n");
    }

  if(strcmp(sfmt,"matroos")==0) {
    fmt=TIMESRIES_FORMAT_MATROOS;
    }

  if(strcmp(sfmt,"classic")==0) {
    fmt=TIMESRIES_FORMAT_CLASSIC;
    }

  if( model == NULL) model=".";

  nst =  mgr_read_stations(input,&sample);

  astro_angles_t astro_angles;
  
  reference=initialize_tide(&astro_angles,start_date);

  if(spectrumFile!=0){
    spectrum = spectrum_init_from_file(spectrumFile,1);
    }
  else if(nbwave==0) {
    i=0;
//  spectrum.waves[i++] = wSa;
//  spectrum.waves[i++] = wSsa;
//  spectrum.waves[i++] = wMm;
//  spectrum.waves[i++] = wMf;

//  spectrum.waves[i++] = wQ1;
//  spectrum.waves[i++] = wO1;
//  spectrum.waves[i++] = wM1;
//    spectrum.waves[nbwave++] = wS1;
//  spectrum.waves[i++] = wPi1;
//  spectrum.waves[i++] = wK1;
//  spectrum.waves[i++] = wKi1;
//  spectrum.waves[i++] = wP1;
  
//  spectrum.waves[i++] = w2N2;
//  spectrum.waves[i++] = wMu2;
//  spectrum.waves[i++] = wN2;
//  spectrum.waves[i++] = wNu2;
//  spectrum.waves[i++] = wM2;
//  spectrum.waves[i++] = wL2;
//  spectrum.waves[i++] = wT2;
    spectrum.waves[nbwave++] = wS2;
//  spectrum.waves[i++] = wK2;
//  spectrum.waves[i++] = wR2;

//  spectrum.waves[i++] = wM4;
    spectrum.n=nbwave;
    }

  separation=printAndGetSpectrumDetails(spectrum);
  printf("Spectrum separation: %g days\n",separation);fflush(stdout);
  if(not isfinite(separation))
    TRAP_ERR_EXIT(-1,"check your spectrum for non separated waves\n");

  if(compare) {
    strcpy(datafile,mgr);
    l=strlen(datafile);
//  work in progress
//     fmt=MGR_FORMAT_IOS;
//     tide_init(&status);
//     mgr_loadstations(datafile, &l, &fmt, &status);
//     mgr_getset(&setptr,&nset,&status);
    }

  if(detiding) {
    a=(double **)malloc(nst*sizeof(double *));
    G=(double **)malloc(nst*sizeof(double *));
    for (k=0;k<nst;k++) {
      a[k]=(double *)malloc(spectrum.n*sizeof(double));
      G[k]=(double *)malloc(spectrum.n*sizeof(double));
      }
    }

  for (k=0;k<nst;k++) {
/* *-----------------------------------------------------------------------------
    read the model elevation file  */
    serie=tseries_t(4);
    if(sample[k].name==0) {
      sprintf(name,"%s/%s%d",model,"sample.",k+1);
      }
    else {
      sprintf(name,"%s/%s%s",model,"sample.",sample[k].name);
      }
    fflush(stdout);
    printf("sea level at %25s: MOD (%s)--",sample[k].name,name);
    status=sample_load(name, serie, &dT[1]);
    if(status ==0) {
      printf("cannot read file for %s, skip station\n",name);
      continue;
      }

/* *----------------------------------------------------------------------------
    limit serie to start and finish date interval*/
//     start=cnes_time(start_date,'d');
//     end  =cnes_time(end_date,'d');

    tseries_t reduced=serie.reduce(start_date,end_date);
    serie.destroy();
    serie=reduced;

    start_date=poctime_getdatecnes(serie.t[0]*24.,'h');
    
    sprintf(name,"%s/%s",path,sample[k].name);
    
    mooring.name=strdup(sample[k].name);
    mooring.lon=sample[k].t;
    mooring.lat=sample[k].p;
    mooring.type=1;
    
    strcpy(tmp,name);
    strcat(tmp,".h");
    status=timeserie_save(tmp, mooring, serie.x[0], serie.mask, serie.t, serie.n, fmt);
//    status=mgr_save_timeserie(string(tmp), mooring, reduced[k][num_step], 'm', format_str, gnu_nice);
    
    strcpy(tmp,name);
    strcat(tmp,".p");
    status=timeserie_save(tmp, mooring, serie.x[3], serie.mask, serie.t, serie.n, fmt);
   
    tseries_t buffer= tseries_t(serie.n, 2);
    valcpy(buffer.t, serie.t, serie.n);
    
/* *----------------------------------------------------------------------------
    process model elevations */
    valcpy(buffer.x[0], serie.x[0], serie.n);
    buffer.mask=serie.mask;
    
    if(test) {
/* *-----------------------------------------------------------------------------
      Create a filtered signal file  */
      strcpy(tmp,name);
      strcat(tmp,".h.1");
      filtering(tmp, buffer.x[0], buffer.t, buffer.mask, buffer.n, dT[1]);
      }
    
    if(doFourier){
      strcpy(tmp,name);
      strcat(tmp,".h.fft");
      fourier(tmp, buffer.x[0], buffer.mask, buffer.n,dT[1], 0);
      }
    
    if(detiding) {
/* *-----------------------------------------------------------------------------
      Detide the sea level signal on request  */
      t0sec=serie.t[0]*24*3600;
      dTsec=dT[1]*24*3600;

      if(sysinit==0) {
        harmonic_init(buffer.x[0], buffer.mask, buffer.n, start_date, dTsec, spectrum, &x);
        sysinit=1;
        }
      if(serie.n!=0)
        harmonic_compute(buffer.x[0], buffer.mask, buffer.x[1], buffer.n, x, spectrum);
//       woce_modelanalysis(buffer,residual,spectrum.wave,&spectrum.n,&t0sec,&dTsec,&nr[1],&status);
      
      status=0;
      for (i=0; i<spectrum.n; i++) {
        zr=buffer.x[0][2*i];
        zi=buffer.x[0][2*i+1];
        a[k][i]=100.*sqrt(zr*zr+zi*zi);
        G[k][i]=atan2(zi,zr)*r2d;
        if(G[k][i]<  0) G[k][i]+=360.0;
        if(G[k][i]>360) G[k][i]-=360.0;
        printf("%4d - %15s, %10s: %6.1f cm %6.1f degrees\n", k,sample[k].name,spectrum.waves[i].name,a[k][i],G[k][i]);
        if(not (a[k][i]<1e5) ){/* can be nan */
          printf("wrong amplitude, most likely because of a poorly separated spectrum\n");
          status++;
          }
        }
      fflush(stdout);
      if(status>0)
        TRAP_ERR_EXIT(-1,"%d wrong amplitudes. See above.\n",status);

      strcpy(tmp,name);
      strcat(tmp,".h");
      status=timeserie_save(tmp, mooring, buffer.x[1], buffer.mask, buffer.t, buffer.n, fmt);
      
      if(test) {
        strcpy(tmp,name);
        strcat(tmp,".h.1");
        filtering(tmp, buffer.x[1], buffer.t, buffer.mask, buffer.n, dT[1]);
        }
      
      if(doFourier){
        strcpy(tmp,name);
        strcat(tmp,".res.fft");
        fourier(tmp, buffer.x[1], buffer.mask, buffer.n, dT[1], 0);
        }
      }

/* *----------------------------------------------------------------------------
    process model ib */
    memcpy(buffer.x[0], serie.x[3], serie.n*sizeof(double));
    
    if(test) {
      strcpy(tmp,name);
      strcat(tmp,".p.1");
      filtering(tmp, buffer.x[0], buffer.t, buffer.mask, buffer.n, dT[1]);
      }
    
    if(doFourier){
      strcpy(tmp,name);
      strcat(tmp,".p.fft");
      fourier(tmp, buffer.x[1], buffer.mask, buffer.n, dT[1], 0);
      }
    
    buffer.destroy();
    serie.destroy();
    }

  if(!compare) TRAP_ERR_EXIT(0,"exiting\n");
/*   out=fopen("score","w"); */

  for (j=0;j<spectrum.n;j++) {
    sprintf(tmp,"%s/score.%s",model,spectrum.waves[j].name);
    out=fopen(tmp,"w");
    for(l=0;l<5;l++) mean[l]=0;
    for(l=0;l<5;l++) rms[l]=0.;
    count=0.;
    for (k=0;k<nst;k++) {
/*  seek the tide gauge index */
      for (i=0;i<nset;i++) {
        n=min(strlen(sample[k].name),8uL);
        if(strncmp(sample[k].name,setptr[i].name,n) == 0) break;
        }
      if(i == nset) {
        break;
        }
      tg=i;
/*  seek the wave index */
      for (i=0;i<setptr[tg].nwave;i++) {
        if(strcmp(spectrum.waves[j].name,setptr[tg].wave[i]) == 0) break;
        }
      if(i == setptr[tg].nwave) {
        break;
        }
      zr=real(setptr[tg].elevation[i]);
      zi=imag(setptr[tg].elevation[i]);
      a1=a[k][j];
      p1=G[k][j];
      a2=100*sqrt(zr*zr+zi*zi);
      p2=360.0-atan2(zi,zr)*r2d;
      if(p2<  0) p2+=360.0;
      if(p2>360) p2-=360.0;
      zr=a[k][j]*cos(-G[k][j]*d2r);
      zi=a[k][j]*sin(-G[k][j]*d2r);
      zr=100* real( setptr[tg].elevation[i] ) -zr;
      zi=100* imag( setptr[tg].elevation[i] ) -zi;
      d=sqrt(zr*zr+zi*zi);
      fprintf(out,"%4d - %6d %25s: %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f %6.1f\n",
            k,setptr[tg].code,setptr[tg].name,a1,a2,a1-a2,p1,p2,p1-p2,d);
      count++;
      l=-1;
      l++;
      d=sqrt((zr*zr+zi*zi)/2.);
      mean[l]=mean[l]+d;
      rms[l]=rms[l]+d*d;
      l++;
      d=a[k][j]-a2;
      mean[l]=mean[l]+d;
      rms[l]=rms[l]+d*d;
      l++;
      d=G[k][j]-p2;
      if(d >  180.) d=d-360;
      if(d < -180.) d=d+360;
      mean[l]=mean[l]+d;
      rms[l]=rms[l]+d*d;
      l++;
      d=zr;
      mean[l]=mean[l]+d;
      rms[l]=rms[l]+d*d;
      l++;
      d=zi;
      mean[l]=mean[l]+d;
      rms[l]=rms[l]+d*d;
      }
    for(l=0;l<5;l++) {
      mean[l]=mean[l]/count;
      rms[l]=sqrt(rms[l]/count - mean[l]*mean[l]);
      }
    fprintf(out,"amplitude bias: %6.1lf, rms: %6.1lf\n",mean[1],rms[1]);
    fprintf(out,"phase bias    : %6.1lf, rms: %6.1lf\n",mean[2],rms[2]);
    fprintf(out,"real part bias: %6.1lf, rms: %6.1lf\n",mean[3],rms[3]);
    fprintf(out,"imag.part bias: %6.1lf, rms: %6.1lf\n",mean[4],rms[4]);
    fprintf(out,"mean deviation: %6.1lf, rms: %6.1lf\n",mean[0],rms[0]);
    fclose(out);
    }

  free(a);
  free(G);
 // free(h);

  TRAP_ERR_EXIT(0,"%s -computation complete ^^^^^^^^^\n",argv[0]);
}
