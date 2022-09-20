#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.h"
#include "archive.h"
#include "poc-time.h"
#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "filter.h"
#include "statistic.h"
#include "functions.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z,velocity[2],elevation,zero=0,beta[3];
  double dt,cut,mask=0;
  double t0,dT,pulsation;
  double sum[2],a1,p1,a2,p2,zr,zi,d;
  float  **h;
  float  *u,*v,*buffer;
  float  count,dummy, hours,start,end,w,dum;
  double *hf,*mhf,*mlf,*lf,*vlf,vhf;
  double time,sampling=0.0,chk01,chk02;
  double *woce_h,*woce_t,woce_mask,dmask;
  int origine,woce_n, nodes[3];
  int fmt,i,j,k,l,L,ndum,rstatus,lfdepth,mfdepth;
  int minstep,maxstep,l1,l2,n,status,nst;
  int  year,month,day,hour;
  int nnde,nstep,elt,m;
  date_t reference;

  size_t size,offset;
  FILE *in,*out,*sout[100];
  FILE *mesh;
  char file1[256],file2[256], wocefile[256], model[256];
  char fout[100][256], outname[256]="\0",rootname[256]="\0";
  pressure_station *sample;
  char *option,*keyword,*s,*meshfile=NULL,*input,*datafile;
  int filtering=0,fix=0,statistic=0;
  int increment=1;
  meta_archive_t info;
  mesh_t benchmark;
  statistic_t stath;

  grid_t *grid;
  pressure_station *setptr;

  fct_echo( argc, argv);

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

        case 's' :
          statistic= 1;
          n++;
          break;

        case 'i' :
          n++;
          sscanf(argv[n],"%d",&increment);
          n++;
          break;


        case 'm' :
          n++;
          meshfile= strdup(argv[n]);
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

  L=strlen(datafile);
  if(strrchr(datafile,(int) '/') != NULL) {
    l=strlen(strrchr(datafile,(int) '/'));
    strncpy(rootname,datafile,L-l+1);
    }
  else  strcpy(rootname,"./");

  printf("treat : %s \n",datafile);

  status=clx_archiveinfo(datafile,&info);
  if(status !=0) goto error;

  if(sampling == 0.) sampling=info.sampling;

  count=info.nframe;
  reference=info.reference;

  status=fe_list(&info.mesh);
/*   if(status !=0) goto error; */


  if(meshfile != NULL) {
/*     fe_freemesh(&info.mesh); */
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&benchmark);
    if(status !=0) goto error;
    status=fe_list(&benchmark);
    if(status !=0) goto error;
    }

  for(m=0; m<info.mesh.ntriangles; m++) {
    status=fe_initaffine(&info.mesh,m);
    }
  nnde=info.mesh.nvtxs;
  for(m=0; m<benchmark.ntriangles; m++) {
    status=fe_initaffine(&benchmark,m);
    }

  if(fix) {
    info.mesh.destroy();
    info.mesh=benchmark;
    }
  for(m=0; m<info.mesh.ntriangles; m++) {
    status=fe_initaffine(&info.mesh,m);
    }

  h=smatrix(0,2,0,nnde-1);
  buffer=(float *) malloc(3*nnde*sizeof(float));

  in=fopen(datafile,"r");
  status=clx_archiveread(in,info,1,h,&time);
/*
  printf("t= %6.1f hours status= %d %s\n",time/3600.0,status, sgetnewdate(info.reference,time));
  fe_integrale01(info.mesh,h[0],sum);
*/
  fclose(in);
  if(fix) status=clx_archivewriteheader(datafile,info,h);

  out=fopen("meanlevel.out","a");

  n=0;
  for (i=0;i<(int) count;i+=increment) {
    in=fopen(datafile,"r");
    status=clx_archiveread(in,info,i+1,h,&time);
    fclose(in);
    if(status != 0) {
      printf("End of file reached... \n");
      break;
      }
/*     if(fmod(i,50) ==0)  */
    fe_integrale01(info.mesh,h[0],sum);
    chk01=sum[0]/sum[1];
    fe_integrale02(info.mesh,h[0],sum,(double) 65.0);
    chk02=sum[0]/sum[1];
    fprintf(out,"%12.1lf %7.4lf %7.4f %s\n",time/24./3600.,chk01,chk02,sgetnewdate(info.reference,time));
    printf("%s :  %7.4lf %7.4lf\n",sgetnewdate(info.reference,time),chk01,chk02);
    if(statistic) {
      stath=get_statistics(h[0],(float) -1000.,nnde);
      stath=get_statistics(h[1],(float) -1000.,nnde);
      stath=get_statistics(h[2],(float) -1000.,nnde); 
      }
    if(fix) {
      for (j=0;j<nnde;j++) {
        h[0][j]=h[0][j]-sum[0]/sum[1];
        }
      status=archive_update2(datafile, info, h, time);
      }
    }

  fclose(out);

  free(h);

/*   fprintf(stderr,"%s -computation complete ^^^^^^^^^\n",argv[0]); */
  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  __ERR_BASE_LINE__("%s -computation aborted ^^^^^^^^^\n",argv[0]);
  exit(-1);
}

