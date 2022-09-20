#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>


#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "archive.h"

#include "filter.h"
#include "statistic.h"

#define nstat 9

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dT;
  double tau;

  int i,j,k,l,nitems;
  int n,status;
  FILE *in;
  char *input=0,*keyword,*pathname=NULL,*s;
  char *onde[10];
  char rootname[256]="\0";
  int nonde;
  date_t start,final,reference,start_date;
  double *time, *z, *residuals, mean;
  double *step, minstep, maxstep;
  statistic_t stat;
  spectrum_t WaveList;
 
  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&start.day,&start.month,&start.year);
          break;
 
       case 'f' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&final.day,&final.month,&final.year);
          break;

        case 'd' :
          i=1;
          while(n+i< argc)  {
            onde[i-1]= strdup(argv[n+i]);
            i++;
            }
          n=n+i;
          nonde=i-1;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        if(input==NULL) {
          input= strdup(argv[n]);
          n++;
          }
        else {
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
          }
        break;
      }
      free(keyword);
    }


  if(pathname==NULL) pathname=strdup(".");
  n= load_serie(input, &time, &z, 'm');

  step=new double[n];
  minstep= 1.e+10;
  maxstep=-1.e+10;

  for (i=0; i<n-1;i++) {
    step[i]=time[i+1]-time[i];
    minstep=MIN(minstep,step[i]);
    maxstep=MAX(maxstep,step[i]);
    }
  stat=get_statistics(step, n-1);

  dT=stat.mean+stat.std/2.0;
  printf("minstep=%lf, maxstep=%lf\n",minstep, maxstep);

/* --------------- init data for harmonic analysis ------------------ */

  nonde=n/3;
  nonde=(time[n-1]-time[0])/dT/2.0;

/* *-----------------------------------------------------------------------
  build a Fourier-like spectrum*/
  nonde=(time[n-1]-time[0])/5.0/2.0;

  WaveList.n   =nonde;
  WaveList.nmax=nonde+10;
  WaveList.waves=new tidal_wave[WaveList.nmax];

  printf("# nbre d'ondes a analyser: %d \n",WaveList.n);

  for (i=0; i<WaveList.n; i++) {
    WaveList.waves[i].omega = (double)(i+1)/(time[n-1]-time[0]);
    sprintf(WaveList.waves[i].name,"%d",i);
    WaveList.waves[i].omega *= 360.0/24.;
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", WaveList.waves[i].name,WaveList.waves[i].omega);
    }

//   nonde=1;
//   WaveList.n=0;
//   WaveList.nmax=nonde+1;
//   WaveList.waves=new tidal_wave[WaveList.nmax];
//   addwave(&WaveList, wM2);
//   for (i=0; i<WaveList.n; i++) {
//     if(WaveList.waves[i].omega == 0.) {
//       WaveList.waves[i].init();
//       }
//     printf ("wave: %10s, pulsation: %12.6f degrees/h \n", WaveList.waves[i].name,WaveList.waves[i].omega);
//     }

/*----------------------------------------------------------------------------*/
//  init_argument(reference);

  residuals=new double[n];

/*-------------------------------------------------------------------------------------
  HARMONIC ANALYSIS */
  status= harmonic_analysis_old(z,time,residuals, &mean, n, WaveList);

/*----free vectors------------------------------------------*/

  goto end;
error:
end:
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
