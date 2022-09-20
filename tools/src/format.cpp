#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "fe.h"
#include "archive.h"
#include "rutin.h"



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *buffer1,*buffer2;
  float **h,*buffer;
  double time,t,previous,sampling,start_time,end_time,start;
  int day,month,year,first_day;
  int i,j,k,l,nndes;
  int n,status,nitems,count,type;
  size_t size;
  FILE *file;
  FILE *in,*out;
  char *datafile,*keyword,*input,*example;
  char *meshfile,*s;
  char *rootname;
  size_t offset;
  date_t date,model_date,reference,start_date,end_date;
  date_t current_date,previous_date;
  meta_archive_t info;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'o' :
          datafile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'e' :
          example= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&reference.day,&reference.month,&reference.year);
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        input= strdup(argv[n]);
        n++;
        break;
      }
      free(keyword);
    }

  status= clx_archiveinfo(example,&info);
  if(status!=0) goto error;
  status=fe_list(&info.mesh);
  if(status!=0) goto error;

  in=fopen(input,"r");
  if(in==NULL) goto error;

  datafile=(char *)malloc(strlen(rootname)+64);

  nndes=info.mesh.nvtxs;

  buffer=(float *)malloc(3*nndes*sizeof(float));
  h=smatrix(0,2,0,nndes-1);

  offset=3*nndes*sizeof(float);

  n=0;
  do {
    size=0;
    if(!(size=fread(&time,sizeof(double),1,in))) {
      printf("End of file reached... \n");
      break;
      }
    if(n%50==0) printf("t= %6.1f hours status= %d %s\n",time/3600.0,status, sgetnewdate(reference,time));
    if(n==0) getnewdate(reference,time,&start_date);
    if(n==1) sampling=time-previous;
    fseek(in,offset,SEEK_CUR);
    previous=time;
    n++;
/*     } while (!feof(in));   */
    } while (n<2);

  info.sampling=sampling;
  info.reference=reference;

  rewind(in);

  size=fread(&time,sizeof(double),1,in);
  size=fread(buffer,sizeof(float),3*nndes,in);

  info.start=time;

  for (j=0;j<nndes;j++) {
    h[0][j]=buffer[3*j];
    h[1][j]=buffer[3*j+1];
    h[2][j]=buffer[3*j+2];
    }

  sprintf(datafile,"%s-%d.%2.2d",rootname,start_date.year,start_date.month);
  status=archive_createheader(datafile,info,h);
  if(status!=0) goto error;

  previous_date=start_date;

  n=1;
  do {
    size=0;
    if(!(size=fread(&time,sizeof(double),1,in))) {
      printf("end of file reached... \n");
      break;
      }
    getnewdate(reference,time,&current_date);
    if(n%10==0)     printf("t= %6.1f hours status= %d %s\n",time/3600.0,status, sgetnewdate(reference,time));
    size=fread(buffer,sizeof(float),3*nndes,in);
    if(size!=3*nndes) goto error;
    for (j=0;j<nndes;j++) {
      h[0][j]=buffer[3*j];
      h[1][j]=buffer[3*j+1];
      h[2][j]=buffer[3*j+2];
      }
    if(current_date.month != previous_date.month) {
      sprintf(datafile,"%s-%d.%2.2d",rootname,current_date.year,current_date.month);
      info.start=time;
      status=archive_createheader(datafile,info,h);
      if(status!=0) goto error;
      }
    else {
      status=archive_update2(datafile, info, h, time);
      if(status!=0) goto error;
      }
    n++;
    previous_date=current_date;
    } while (!feof(in));


  free(buffer);
  free_smatrix(h,0,2,0,nndes-1);
  fclose(in);
  fclose(out);
  printf("end of format ... \n");


 error:
  free(buffer);
  free_smatrix(h,0,2,0,nndes-1);
  fclose(in);
  fclose(out);
  printf("format aborted ... \n");
}
