#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "archive.h"
#include "poc-time.h"
#include "functions.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float **buffer1,**buffer2,**buffer3;

  float  count, hours,start,end,ibd;
  double time1,time2,time3,restart_time,dt;
  int day,month,year;
  int nndes,i,j,k,l,L,n,status,nread;
  size_t size;
  FILE *in,*out;
  char *hfile,*keyword,*example;
  char *s,*sdate;
  char filename[256],line[256];
  char rootname[256]="\0",tmp[256];
  date_t model_date,date,reference;
  meta_archive_t info;
  state2D_t *set;

  fct_echo(argc,argv);
  printf(" ^^^^^^^^^^^^^ start computation ^^^^^^^^^^^^^\n");

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'h' :
          hfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'e' :
          example= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&day,&month,&year);
          break;

        case 'i' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%lf",&dt);
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

  L=strlen(hfile);
  if(strrchr(hfile,(int) '/') != NULL) {
    l=strlen(strrchr(hfile,(int) '/'));
    strncpy(rootname,hfile,L-l+1);
    }
  else  strcpy(rootname,"./");

  printf("rootname: %s \n",rootname);

  status=clx_archiveinfo(hfile,&info);

  reference=info.reference;
  nndes=info.mesh.nvtxs;

  restart_time=((mjd(month,day,year)-mjd(reference.month,reference.day,reference.year))*24.)*3600.;
  printf("restart= %6.1f hours status= %d %s\n",restart_time/3600.0,status, sgetnewdate(info.reference,restart_time));
  printf("dt= %6.1f deconds\n",dt);

  in=fopen(hfile,"r");

  buffer1=smatrix(0,2,0,nndes-1);
  buffer2=smatrix(0,2,0,nndes-1);
  buffer3=smatrix(0,2,0,nndes-1);

  i=(int) ((restart_time-info.start)/info.sampling);

  status=clx_archiveread(in,info,i,buffer1,&time1);
  printf("t= %6.1f hours status= %d %s\n",time1/3600.0,status, sgetnewdate(info.reference,time1));
  status=clx_archiveread(in,info,i+1,buffer2,&time2);
  printf("t= %6.1f hours status= %d %s\n",time2/3600.0,status, sgetnewdate(info.reference,time2));
  status=clx_archiveread(in,info,i+2,buffer3,&time3);
  printf("t= %6.1f hours status= %d %s\n",time3/3600.0,status, sgetnewdate(info.reference,time3));

  fclose(in);

  set= (state2D_t *) malloc(nndes*sizeof(state2D_t));

  in=fopen(example,"r");

  fgets(line,256,in);
  for(i = 1;i<=nndes;i++) {
    fgets(line,256,in);
    if(feof(in)) goto error;
    nread=sscanf(line," %lf %lf %lf %f %f %lf %lf %lf %lf %lf\n", &ST(set,h,n),
                    &ST(set,H,i), &ST(set,Hm,i),
                    &ST(set,u,i), &ST(set,v,i),
                    &ST(set,hmean,i),
                    &ST(set,umean,i),
                    &ST(set,vmean,i),
                    &ST(set,humean,i),
                    &ST(set,hvmean,i));
    if(nread<5) goto error;
    }

  fclose(in);

  for (i=0;i<nndes;i++) {
    n=i+1;
    ST(set,h,n) =info.mesh.vertices[i].h;
    ST(set,H,n) =buffer2[0][i]+ST(set,h,n);
    ST(set,Hm,n)=ST(set,H,n)-0.5*(dt/info.sampling)*(buffer2[0][i]-buffer1[0][i])
                            -0.5*(dt/info.sampling)*(buffer3[0][i]-buffer2[0][i]);

/*     a=0.5*(buffer1[0][i]+buffer3[0][i])-buffer2[0][i]; */
/*     b=0.5*(buffer3[0][i]-buffer1[0][i]); */

/*     ST(set,Hm,n)=ST(set,H,n)-(2*a*(-dt/info.sampling)+b)*dt/info.sampling; */
/*     ST(set,Hm,n)=ST(set,H,n)-(a*(-dt/info.sampling)+b)*dt/info.sampling; */

    ST(set,u,n) =buffer2[1][i];
    ST(set,v,n) =buffer2[2][i];
    }

/*----------------------------------------------------------------------
  save model state for further restart*/

  sdate=sgetnewdate(info.reference,time2);
  getnewdate(info.reference,time2,&date);
  printf("#Continuation in %s at %s (t= %f hrs)\n","restart",sdate,time2/3600.0);


  sprintf(filename,"restart.%2.2d-%2.2d-%4.4d",date.day,date.month,date.year);

  if ((out = fopen(filename, "w")) == NULL) {
    __ERR_BASE_LINE__("exiting\n");exit(-1);
    }

  fprintf(out,"%f (model time), h,u,v at %s (t= %f hrs)\n",time2,sdate,time2/3600.);
  for(i=1;i<=nndes;i++) {
    fprintf(out," %lf %lf %lf %f %f  %lf %lf %lf %lf %lf\n", ST(set,h,i),
                    ST(set,H,i), ST(set,Hm,i),
                    ST(set,u,i), ST(set,v,i),
                    ST(set,hmean,i),
                    ST(set,umean,i),
                    ST(set,vmean,i),
                    ST(set,humean,i),
                    ST(set,hvmean,i));
    }
  fclose(out);

  free(sdate);

  free_smatrix(buffer1,0,2,0,nndes-1);
  free_smatrix(buffer2,0,2,0,nndes-1);

  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  gmerror("read error");
  __ERR_BASE_LINE__("exiting\n");exit(-1);
  
}
