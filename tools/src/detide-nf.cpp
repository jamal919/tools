#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "tools-structures.h"
     
#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "archive.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z,dum[3],velocity[2],elevation;
  double t0,dT,pulsation;
  double *serie[500],a1,p1,a2,p2,zr,zi,d;

  float  **h;
  float  *u,*v,count,dummy, hours,start;
  float  mean,rms;
  double time;
  int fmt,i,j,k,l,ndum,rstatus;
  int minstep,maxstep,l1,l2,n,status,nst;
  size_t size;
  FILE *in,*out;
  char *option,*keyword,*s,*meshfile=NULL,*output,*datafile,*mgrfile;
  int filtering=0,fix=0;
  meta_archive_t info;

  //char spectrum_t[10]={"S1       "};
  int  nwave=1;

  float dtr=M_PI/180.0;

  grid_t *grid;
  pressure_station *setptr;

  fct_echo(argc,argv);
 
  fprintf(stderr,"%s  -starting computation *********\n",argv[0]);
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
        case 'f' :
          fix= 1;
          n++;
          break;

        case 'm' :
          n++;
          meshfile= strdup(argv[n]);
          n++;
          break;

        case 'o' :
          n++;
          output= strdup(argv[n]);
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        datafile= strdup(argv[n]);
        n++;
        break;
      }
      free(keyword);
    }


  status=clx_archiveinfo(datafile,&info);
  if(status !=0) goto error;

  dT=info.sampling;

  count=info.nframe;
  
  mgrfile=strdup("/data/bartok2/DataBase/Tidal_Gauges/mgr/westmed.mgr");
  
  l=strlen(datafile);

  fmt=0;
  /* to be replaced with C routines
  tide_init(&status);

  mgr_loadstations(mgrfile, &l, &fmt, &rstatus);
  mgr_getset(&setptr,&nst,&status);
  */
  for (k=0;k<nst;k++)serie[k]= (double *) malloc( (long unsigned int)(count) *sizeof(double));
  
  h=smatrix(0,2,0,info.mesh.nvtxs-1);
  n=0;
  for (i=1;i<=count;i++) {
    in=fopen(datafile,"r");
    status=clx_archiveread(in,info,i,h,&time);
    fclose(in);
    if(status != 0) {
      printf("End of file reached... \n");
      break;
      }
    if(time > 3*24*3600.) {
      if(n == 0) t0=time;
      n++;
      for (k=0;k<nst;k++)
        serie[k][n-1] = h[0][setptr[k].code-1];
      }
    }
  
  out=fopen(argv[2],"w");

  mean=0;
  rms=0.;

  for (k=0;k<nst;k++) {
  /* to be replaced with C routines
    woce_modelanalysis(serie[k],spectrum_t,&nwave,&t0,&dT,&n,&rstatus);
  */
    zr=serie[k][0];
    zi=serie[k][1];
    a1=100.0*sqrt(zr*zr+zi*zi);
    p1=360.0-atan2(zi,zr)/dtr;
    zr=real(setptr[k].elevation[0]);
    zi=imag(setptr[k].elevation[0]);
    a2=100*sqrt(zr*zr+zi*zi);
    p2=360.0-atan2(zi,zr)/dtr;
    zr=100*( real(setptr[k].elevation[0]) -serie[k][0]);
    zi=100*( imag(setptr[k].elevation[0]) -serie[k][1]);
    d=sqrt(zr*zr+zi*zi);
    fprintf(out,"%4d - %6d %25s: %6.1f %6.1f %6.1f %6.1f %6.1f\n",
            k,setptr[k].code,setptr[k].name,a1,a2,p1,p2,d);
    d=d;
    mean=mean+d;
    rms=rms+d*d;
   }
  for (k=0;k<nst;k++) free(serie[k]);

  mean=mean/nst;
  rms=sqrt(rms/nst - mean*mean);
  fprintf(out,"mean deviation: %6.1f, rms: %6.1f\n",mean,rms);

  free(h);
  
  __ERR_BASE_LINE__("%s -computation complete ^^^^^^^^^\n",argv[0]);
  exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
  
}
