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
#include "poc-time.h"
#include "mgr.h"
#include "functions.h"


#define nstat 9

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double wave_clock(vector<mgr_t> mgr[2],spectrum_t WaveList, char wave[2][10])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  h = a * cos(omega*t+V0-G)

  V = omega*t+V0

  In case of initial time error:

  true(V0-G)=false(V0-G)

  dV0=true(V0)-false(V0) + n * 2pi = omega * (true(t)-false(t))
                                   = omega * dt
                                   = true(G)-false(G) + n * 2pi
                                   = dG' + n * 2pi

  omega*dt=dG'+ n * 2 * pi = omega*dt' + n * 2 * pi

  dt=dt'+ n * T

  approximation: find shift within 15 days lag

  semi-diurnal :  Ns
  diurnal      :  Nd

  01,K1 : 1° per hour
  M2,S2 : 1° per hour

  First, estimate N:

  Ns (T(M2)-T(S2)) = -dt'(M2)-dt'(S2)
  Nd (T(K1)-T(O1)) = -dt'(K1)-dt'(O1)

  Second, seek for best fit:

  dt=average(dt'+ n * T)

-----------------------------------------------------------------------------*/
{
int  k[2][2];

  int i,j,w;

  i=0;
  k[0][0]=mgr[0][i].wave_index(wave[0]);
  k[1][0]=mgr[0][i].wave_index(wave[1]);
  k[0][1]=mgr[1][i].wave_index(wave[0]);
  k[1][1]=mgr[1][i].wave_index(wave[1]);

  i=0;
  for (j=0;j<mgr[0][i].nwave;j++) {
    WaveList.waves[j]=wave_from_code(mgr[0][i].data[j].constituent.code);
    }

double dG[2];
double dtt,dt[2],T[2];

  for(w=0;w<2;w++) {
    WaveList.waves[k[w][0]].init();
/* *----------------------------------------------------------------------------
    compute phase lag difference*/
    dG[w]=mgr[1][i].data[k[w][1]].phi-mgr[0][i].data[k[w][0]].phi;
    T[w]=360./WaveList.waves[k[w][0]].omega;
    dt[w]=dG[w]/WaveList.waves[k[w][0]].omega;
    }

int N;

  for(N=-5;N<5;N++) {
    dtt=0.5*(dt[0]+N*T[0]+dt[1]+N*T[1]);
    printf("%s --- %s %s\n",fct_hms(dtt),fct_hms(dt[0]+N*T[0]),fct_hms(dt[1]+N*T[1]));
    }

  N= -floor((dt[0]-dt[1])/(T[0]-T[1])+0.5);

  printf("%best fit---------------------------------------\n");
  dtt=0.5*(dt[0]+N*T[0]+dt[1]+N*T[1]);
  printf("%s --- %s: %s %s: %s\n",fct_hms(dtt),wave[0],fct_hms(dt[0]+N*T[0]),wave[1],fct_hms(dt[1]+N*T[1]));
  printf("%best fit---------------------------------------\n");

  return(dtt);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int clock(vector<mgr_t> mgr[2],spectrum_t WaveList)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

double dt;
char wave[2][10];

  strcpy(wave[0],"M2");
  strcpy(wave[1],"S2");

  dt=wave_clock(mgr,WaveList,wave);

  strcpy(wave[0],"M2");
  strcpy(wave[1],"N2");

  dt=wave_clock(mgr,WaveList,wave);

  strcpy(wave[0],"S2");
  strcpy(wave[1],"N2");

  dt=wave_clock(mgr,WaveList,wave);

  strcpy(wave[0],"2N2");
  strcpy(wave[1],"K2");

  dt=wave_clock(mgr,WaveList,wave);

  strcpy(wave[0],"K1");
  strcpy(wave[1],"O1");

  dt=wave_clock(mgr,WaveList,wave);
}




/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dT;
  double tau;

  int i,j,k,l,nitems;
  int n,status;
  FILE *in;
  char *input=0,*keyword,*pathname=NULL,*reference=NULL,*s;
  char *onde[10];
  char rootname[256]="\0";
  int nonde;
  date_t start,final,start_date;
  double dt;
  statistic_t stat;
  spectrum_t WaveList;

  vector<mgr_t> mgr[2];
  int nmgr[2];

  fct_echo( argc, argv);

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

        case 'r' :
          reference= strdup(argv[n+1]);
          n++;
          n++;
          break;

       case 'f' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&final.day,&final.month,&final.year);
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

  WaveList.waves=new tidal_wave[100];
  WaveList.nmax=100;

  nmgr[0]=mgr_load(input,      mgr[0]);
  nmgr[1]=mgr_load(reference,  mgr[1]);


  dt=clock(mgr,WaveList);

  goto end;
error:
end:
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
