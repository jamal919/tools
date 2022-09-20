

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
#include "functions.h"
#include "matrix.h"

#define nstat 9

/*----------------------------------------------------------------------

Compute global statistics for the MOG2D outputs

----------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void smagorP1NQ(state2d_obsolete_t * P1_state, parameter_obsolete_t * P1_data, int nndes,
                  triangle_t * elt, int nelt, double Ah)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------------*/

  int i, j, k, l, m, n;
  int n1, n2, n3;
  double D1, D2, AH, mean, rms, A, factor;
  double ux, uy, vx, vy, tmp[nelt];

/*----------------------------------------------------------------------
  FACE VALUE (QUODDY 4 ???)
  double SMAGR  =   0.28;

  d/dx= DQ/2A (m^-1)
  d/dy=-DP/2A

  D1 and D2 are s^-1
  AH and Asmag are m^2/s
----------------------------------------------------------------------*/

  double SMAGR = 0.28;
//  double new,damp;
  double damp;

  for(n = 0; n < nndes; n++) {
    P1_state[n].Ah = 0.0;
    P1_state[n].ux = 0.0;
    P1_state[n].uy = 0.0;
    P1_state[n].vx = 0.0;
    P1_state[n].vy = 0.0;
    }

#pragma omp parallel
#pragma shared (elt,P1_state,P1_data,nelt,SMAGR,smb1)
#pragma local (l,n1,n2,n3,AH)
#pragma local (D1,D2,A)

/*BRACKET FOR PARALLEL COMPUTATION*/
  {
#pragma pfor iterate (l=1;nelt;1)
/*----------------------------------------------------------------------
    Use of variational approach in the nodal quadrature frame*/
    for(l = 0; l < nelt; l++) {
      A = elt[l].Area;
      n1 = elt[l].vertex[0];
      n2 = elt[l].vertex[1];
      n3 = elt[l].vertex[2];
      ux =  (P1_state[n1].u * elt[l].DQ[0] / P1_data[n1].C
           + P1_state[n2].u * elt[l].DQ[1] / P1_data[n2].C
           + P1_state[n3].u * elt[l].DQ[2] / P1_data[n3].C) / (2. * A);
      uy = -(P1_state[n1].u * elt[l].DP[0]
           + P1_state[n2].u * elt[l].DP[1]
           + P1_state[n3].u * elt[l].DP[2]) / (2. * A);
      vx =  (P1_state[n1].v * elt[l].DQ[0] / P1_data[n1].C
           + P1_state[n2].v * elt[l].DQ[1] / P1_data[n2].C
           + P1_state[n3].v * elt[l].DQ[2] / P1_data[n3].C) / (2. * A);
      vy = -(P1_state[n1].v * elt[l].DP[0]
           + P1_state[n2].v * elt[l].DP[1]
           + P1_state[n3].v * elt[l].DP[2]) / (2. * A);

      for(j = 0; j < 3; j++) {
        n = elt[l].vertex[j];
        P1_state[n].ux += ux * A / 3;
        P1_state[n].uy += uy * A / 3;
        P1_state[n].vx += vx * A / 3;
        P1_state[n].vy += vy * A / 3;
        }
/*----------------------------------------------------------------------
      grad u/x -grad v/y*/
      D1 = ux - vy;
/*----------------------------------------------------------------------
      grad v/x +grad u/y*/
      D2 = vx + uy;
      AH = SMAGR * 2.0 * A * sqrt(D1 * D1 + D2 * D2);
      tmp[l] = AH * A / 3;
      }
#pragma synchronize
/*BRACKET FOR PARALLEL COMPUTATION*/
    }

  for(l = 0; l < nelt; l++) {
    n1 = elt[l].vertex[0];
    n2 = elt[l].vertex[1];
    n3 = elt[l].vertex[2];
    P1_state[n1].Ah += tmp[l];
    P1_state[n2].Ah += tmp[l];
    P1_state[n3].Ah += tmp[l];
    }

  for(n = 0; n < nndes; n++) {
    factor=1./P1_data[n].lmm;
    P1_state[n].ux *= factor;
    P1_state[n].uy *= factor;
    P1_state[n].vx *= factor;
    P1_state[n].vy *= factor;
    P1_state[n].Ah *= factor;
    if(P1_state[n].Ah < Ah)
      P1_state[n].Ah = Ah;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void WE_viscosityLGP1(state2d_obsolete_t * P1_state, parameter_obsolete_t * P1_data,
                            int nnde, triangle_t * elt, int nelt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*----------------------------------------------------------------------------*/
  int i, j, l, n, k, m;
  int cnt, *list;
  double CC, CC1, uDQ, uDP, vDQ, vDP;
  double u, v, A;
/*   triangle_t element; */
/*----------------------------------------------------------------------------*/

  for(n = 0; n < nnde; n++) {
    P1_state[n].Ex = 0.0;
    P1_state[n].Ey = 0.0;
    }

  for(l = 0; l < nelt; l++) {
    CC = 0.0;
    CC1 = 0.0;
    uDQ = 0.0;
    uDP = 0.0;
    vDQ = 0.0;
    vDP = 0.0;
    A = elt[l].Area * (double) 12.;
    for(j = 0; j < 3; j++) {
      n = elt[l].vertex[j];
      CC += P1_data[n].C;
      CC1 += 1. / P1_data[n].C;
      u = (double) P1_state[n].u;
      v = (double) P1_state[n].v;
      uDQ += u * elt[l].DQ[j];  /* = 2*A dudt */
      uDP += u * elt[l].DP[j];
      vDQ += v * elt[l].DQ[j];
      vDP += v * elt[l].DP[j];
      }
    for(j = 0; j < 3; j++) {
      n = elt[l].vertex[j];
      P1_state[n].Ex +=
            (float) (elt[l].DQ[j] * uDQ * CC1 + elt[l].DP[j] * uDP * CC) / A;
      P1_state[n].Ey +=
            (float) (elt[l].DQ[j] * vDQ * CC1 + elt[l].DP[j] * vDP * CC) / A;
      }
    }

/*----------------------------------------------------------------------
  lmm*C= sum over elts of (AC/3)*/
  for(n = 0; n < nnde; n++) {
    P1_state[n].Ex *= -P1_state[n].Ah / (P1_data[n].lmm * P1_data[n].C);
    P1_state[n].Ey *= -P1_state[n].Ah / (P1_data[n].lmm * P1_data[n].C);
  }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z,dum[3],dum1[3],dum2[3],velocity[2],elevation,ws[2],pressure;
  float *buffer[10],*a,*G;
  double origine,dT;
  double tau;

  float *h,*hmean;
  double time;
  int count;
  int fmt,i,j,k,l,L,ndum,rstatus,shift,nitems;
  int n,status;
  int day,month,year,nday,loop,elapsed,step;
  FILE *in;
  char hfile[256],pfile[256],*hroot,*proot,*keyword,*zone,*s,*choice,*pathname=NULL;
  char *onde[11],file1[256],file2[256],out[256];
  char rootname[256]="\0",tmp[256],output[256];
  grid_t grid;
  int nonde;
  date_t start,final,reference,start_date;
  meta_archive_t info;
  state2d_obsolete_t   *state;
  parameter_obsolete_t *data;

  int nndes;
  spectrum_t WaveList, AnalysisList;
  date_t lastref;

  fct_echo( argc, argv);

  i=0;
  onde[i++]=strdup("Z0");
  onde[i]=NULL;

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {

        case 'p' :
          pathname= strdup(argv[n+1]);  /* directory */
          n++;
          n++;
          break;

        case 'h' :
          hroot= strdup(argv[n+1]);    /* (analysis) */
          n++;
          n++;
          break;

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
          n++;
          i=pos((char*)NULL,onde,11);
          for(;n<argc;i++,n++) {
            onde[i]= strdup(argv[n]);
            }
          onde[i]=NULL;
          nonde=i;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        __OUT_BASE_LINE__("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
      free(keyword);
    }


  if(pathname==NULL) pathname=strdup(".");

/*---------------- read the mesh file (neigbour format) --------------*/
/*  sprintf(hfile,"%s-%4.4d.%2.2d.reduced-6",hroot,start.year,start.month);*/
  sprintf(hfile,"%s/%s-%4.4d.%2.2d",pathname,hroot,start.year,start.month);

  printf("get info from file %s\n",hfile);

  status= clx_archiveinfo(hfile,&info);
  if(status !=0) {
    __OUT_BASE_LINE__("unable to read properly file %s; abort ...\n",hfile);
    exit(-1);
    }
  status=fe_list(&info.mesh);

  nndes=info.mesh.nvtxs;
  printf("nndes=%d\n",nndes);

  dT=info.sampling; /* dT en secondes     */
  reference=info.reference;


 /*-------------------- allocate memory for vectors -------------------*/
  for(k=0;k<5;k++) {
    buffer[k]=(float *) malloc(nndes*sizeof(float));
    if(buffer[k] == NULL) {
      printf("#memory allocation error for buffer N= %d \n",nndes);
      goto error;
      }
    }

  h=(float *) malloc(nndes*sizeof(float));
  if(h == NULL) goto error;
  hmean=(float *) malloc(nndes*sizeof(float));
  if(hmean == NULL) goto error;

  state=new state2d_obsolete_t[nndes];
  data =new parameter_obsolete_t[nndes];

/* --------------- init data for harmonic analysis ------------------ */
/*ATTENTION: l'egalite de structures ne marche pas bien sur pccls ... */

  WaveList.init(initialize_tide(),onde);
  printf("# nbre d'ondes a analyser: %d \n",WaveList.n);

  for (i=0; i<WaveList.n; i++) {
    if(WaveList.waves[i].omega == 0.) {
      WaveList.waves[i].init();
      }
    printf ("wave: %10s, pulsation: %12.6f degrees/h \n", WaveList.waves[i].name,WaveList.waves[i].omega);
    }

  for (i=0; i<WaveList.n; i++) {
    for (j=i+1; j<WaveList.n; j++) {
      tau=deltaOmega2separation(WaveList.waves[j].omega-WaveList.waves[i].omega);
      printf ("wave: %10s %10s, separation: %9.3f days \n",
        WaveList.waves[i].name,WaveList.waves[j].name,tau);
      }
    }

 /*-----------------------------------------------------------------------*/
  astro_angles_t astro_angles;
  init_argument(&astro_angles,reference);

  for(i=0;i<nndes;i++) {
    hmean[i]=0.;
    }

  archive_freeinfo(&info);

  harmonic_init(WaveList,AnalysisList,nndes);
  count=0;

  for (year=start.year;year<=final.year;year++) {
    for (month=1;month<=12;month++) {
      if ((year==start.year) && (month<start.month)) continue;
      if ((year==final.year) && (month>final.month)) break ;
      sprintf(hfile,"%s/%s-%4.4d.%2.2d",pathname,hroot,year,month);
      printf("read %s  ...\n",hfile);
/*-------------------------------------------------------------------------------------
      read the output files (MOG2D format)*/

      status= clx_archiveinfo(hfile,&info);
      if(status !=0) goto read_error;
      printf("info.nframe= %d\n",info.nframe);

      in=fopen(hfile,"r");
      for(step=1;step<=info.nframe;step++) {
/*         if (count==48) break ; */
        status=clx_archiveread(in, info, step, buffer, &time);
        if(status !=0) goto read_error;
        if(step%25 ==0) printf("t= %6.1f hours frame= %d %s\n",time/3600.0,step,sgetnewdate(info.reference,time));
        for (n=0;n<nndes;n++) {
          state[n].H=buffer[0][n];
          state[n].u=buffer[1][n];
          state[n].v=buffer[2][n];
          }
/*
        smagorP1NQ(state,data,info.mesh.nvtxs,info.mesh.triangles,info.mesh.ntriangles,(double) 1.0);
        WE_viscosityLGP1(state,data,info.mesh.nvtxs,info.mesh.triangles,info.mesh.ntriangles);
        for (n=0;n<nndes;n++) {
          buffer[3][n]=state[n].Ex;
          buffer[4][n]=state[n].Ey;
          }
*/
        harmonic_storage_obsolete(AnalysisList, time,nndes,hmean,buffer,&count, astro_angles);
        }
      fclose(in);
      archive_freeinfo(&info);
      }
    }

/*-------------------------------------------------------------------------------------
  HARMONIC ANALYSIS */

  harmonic_end(AnalysisList, start,final,nndes,count);

/*----free vectors------------------------------------------*/
   for(k=0;k<3;k++) {
    free(buffer[k]);
    }

  free(h);
  free(hmean);

  goto end;

error:
  printf(" error in allocating memory  ...\n");
  read_error:
  printf(" error in reading %s ...\n",hfile);
  free(h);
  free(hmean);
  __ERR_BASE_LINE__("exiting\n");exit(-1);

end:
  __OUT_BASE_LINE__(" ^^^^^^^^^^^^^ end of computation ^^^^^^^^^^^^^\n");
  exit(0);
}
