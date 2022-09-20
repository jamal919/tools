#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
     
#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "fe.h"
#include "poc-time.h"
#include "archive.h"
#include "bmg.h"
#include "map.h"

#define nstat 9


/*----------------------------------------------------------------------

Compute global statistics for the MOG2D outputs

----------------------------------------------------------------------*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_zonegrid_fatal(char *zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;
  int zone_initialised=0;

/* NEA grid */
  if(strcmp(zone,"NEA-LR")==0) {
  map_set2Dgrid(&grid,-20.0,+29.5,+15.0,+65.0,0.5,0.5);
  grid.nx  = 71;
  grid.ny  = 72;
  zone_initialised=1;
  }
  
  if(strcmp(zone,"NEA")==0) {
  map_set2Dgrid(&grid,-20.0,+29.5,+15.0,+65.0,0.25,0.25);
  grid.nx  = 141;
  grid.ny  = 143;
  zone_initialised=1;
  }
  
  if(strcmp(zone,"NEA-HR")==0) {
  map_set2Dgrid(&grid,-20.0,+29.5,+15.0,+65.0,0.1,0.1);
  grid.nx  = 351;
  grid.ny  = 356;
  zone_initialised=1;
  }
  
/* kerguelen grid */
  if(strcmp(zone,"kerguelen")==0) {
  grid.xmin= +40.0 ;
  grid.ymin= -57.5;
  grid.dx  =   0.5;
  grid.dy  =   0.5;
  grid.nx  = 96;
  grid.ny  = 51;
  zone_initialised=1;
  }
  
/* indian ocean grid */
  if(strcmp(zone,"indian")==0) {
  grid.xmin= +10.0 ;
  grid.ymin= -75.0;
  grid.dx  =   1.0;
  grid.dy  =   1.0;
  grid.nx  = 141;
  grid.ny  = 111;
  zone_initialised=1;
  }
  
/* global grid */
  if(strcmp(zone,"global")==0) {
  grid.xmin=0.;
  grid.ymin=-90.;
  grid.dx=1;
  grid.dy=1;
  grid.nx=361;
  grid.ny=181;
  zone_initialised=1;
  }
  
/* global grid for R. ponte */
  if(strcmp(zone,"global-RP")==0) {
  grid.xmin=0.;
  grid.ymin=-90.;
  grid.dx=1.125;
  grid.dy=1.125;
  grid.nx=320;
  grid.ny=161;
  zone_initialised=1;
  }
  
/* mesdsea grid, low resolution */
  if(strcmp(zone,"medsea-LR")==0) {
  grid.xmin= -10.0 ;
  grid.ymin=  27.5;
  grid.dx  =   0.5;
  grid.dy  =   0.5;
  grid.nx  = 101;
  grid.ny  = 41;
  zone_initialised=1;
  }
  
/* mesdsea grid, high resolution */
  if(strcmp(zone,"medsea-HR")==0) {
  map_set2Dgrid(&grid,-10.0,27.5,40.0,47.5,0.1,0.1);
  grid.nx  = 501;
  grid.ny  = 201;
  zone_initialised=1;
  }
  
/* mesdsea grid, normal resolution */
  if(strcmp(zone,"medsea")==0) {
  map_set2Dgrid(&grid,-10.0,27.5,40.0,47.5,0.25,0.25);
  grid.nx  = 201;
  grid.ny  = 81;
  zone_initialised=1;
  }
  if(zone_initialised!=1) {__ERR_BASE_LINE__("exiting\n");exit(-1);}

  return(grid);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z,dum[3],dum1[3],dum2[3],velocity[2],elevation,ws[2],pressure;
  float **buffer1,**buffer2,**buffer3,**buffer4,*buf,*bufx,*bufy;
  float **swap;

  float  H,h,u,v,p,hu,hv,wsx,wsy,hwsx,hwsy,*avg[nstat],*rms[nstat];
  float  count, hours,mask=1.e+10,ibd;
  double time,t,t0,model_start,decal=0.;
  int step,nndes,incr,decald;
  int d1,d2,i,j,k,l,L,init=0;
  int l1,l2,n,status;
  int day,month,year,nday,loop,elapsed;
  int *element;
  bmgheader_t header;

  size_t size;
  FILE *in1,*in2,*out;
  char hfile[1024],pfile[1024];
  char *hroot=NULL,*proot=NULL,*keyword,*zone,*s;

  char file1[256],file2[256];
  char rootname[256]="\0",tmp[256],output[1024],input[1024];
  grid_t grid;
  meta_archive_t info;
  date_t actual;

  date_t start,final;

  fct_echo(argc,argv);
 
  n=1;
  while (n < argc)
    {
    keyword=strdup(argv[n]);
    switch (keyword[0])
      {
      case '-':
        switch (keyword[1])
        {
        case 'p' :
          proot= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'h' :
          hroot= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'l' :
          s= strdup(argv[n+1]);
          sscanf(s,"%d",&nday);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;


        case 's' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&start.day,&start.month,&start.year);
          break;

        case 'e' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&final.day,&final.month,&final.year);
          break;

        case 'd' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%lf",&decal);
          break;
/*
        case 'l' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&day,&month,&year);
          break;
*/

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

  if (strrchr(hfile,(int) '/') != NULL) {
    l=strlen(strrchr(hfile,(int) '/'));
    L=strlen(hfile);
    strncpy(rootname,hfile,L-l+1);
    }
  else strcpy(rootname,"./");

  printf("rootname: %s \n",rootname);

  decald=(int)(decal+0.25);
  printf("decal day=%lf (%d)\n",decal,decald);
  
 /*---------------- read the mesh file (neigbour format) --------------*/

  if(hroot==NULL) goto error;
  if(proot==NULL) goto error;

  sprintf(hfile,"%s-%4.4d.%2.2d",hroot,start.year,start.month);
  status= clx_archiveinfo(hfile,&info);
  if(status!=0) goto error;
  status=fe_list(&info.mesh);
  if(status!=0) goto error;

  nndes=info.mesh.nvtxs;

 /*-------------------- allocate memory for vectors -------------------*/

  for(k=0;k<nstat;k++) avg[k]=(float *)malloc(nndes*sizeof(float));
  for(k=0;k<nstat;k++) rms[k]=(float *)malloc(nndes*sizeof(float));
 
  buffer1=(float **) malloc(3*sizeof(float *));
  buffer2=(float **) malloc(3*sizeof(float *));
  buffer3=(float **) malloc(3*sizeof(float *));
  buffer4=(float **) malloc(3*sizeof(float *));

  for(k=0;k<3;k++) {
    buffer1[k]=(float *)malloc(nndes*sizeof(float));
    if(buffer1[k] == NULL) {
      printf("#memory allocation error for buffer N= %d \n",nndes);
      goto error;
      }
    buffer2[k]=(float *)malloc(nndes*sizeof(float));
    if(buffer2[k] == NULL) {
      printf("#memory allocation error for buffer N= %d \n",nndes);
      goto error;
      }
    buffer3[k]=(float *)malloc(nndes*sizeof(float));
    if(buffer3[k] == NULL) {
      printf("#memory allocation error for buffer N= %d \n",nndes);
      goto error;
      }
    buffer4[k]=(float *)malloc(nndes*sizeof(float));
    if(buffer3[k] == NULL) {
      printf("#memory allocation error for buffer N= %d \n",nndes);
      goto error;
      }
    }

  for(loop=0;loop<1;loop++) {
    for (i=0;i<nndes;i++) {
      for(k=0;k<nstat;k++) avg[k][i]=0.;
      for(k=0;k<nstat;k++) rms[k][i]=0.;
      };

  grid=get_zonegrid_fatal(zone);
  status=map_completegridaxis_2(&grid);

  element=fe_scan_elements(info.mesh,grid,0);
  if (element == NULL) {__ERR_BASE_LINE__("exiting\n");exit(1);}

  buf=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  if (buf == NULL) {__ERR_BASE_LINE__("exiting\n");exit(2);}
  printf("%4.4d.%2.2d\n",start.year,start.month);
  printf("%4.4d.%2.2d\n",final.year,final.month);

  d1=julian_day(start.month, 1,start.year)-julian_day(1,1,1950);
  d2=julian_day(final.month,30,final.year)-julian_day(1,1,1950);

  count=0;
  for (day=d1;day<=d2-decald;day+=10) {
    calendary(day+decald,&actual);

    sprintf(hfile,"%s-%4.4d.%2.2d",hroot,actual.year,actual.month);
    sprintf(pfile,"%s-%4.4d.%2.2d",proot,actual.year,actual.month);

    in1=fopen(hfile,"r");
    in2=fopen(pfile,"r");
  
    printf("read %s & %s ...\n",hfile,pfile);

 /*---------------- read the output files (MOG2D format) --------------*/

    status= clx_archiveinfo(hfile,&info);

    t0=(julian_day(info.reference.month,info.reference.day,info.reference.year)-CNES0jd)*24.*3600.+info.reference.second;

    model_start=t0+info.start;
    if(decal != 0.)
      t=floor((day+decal)*24.*3600.+0.5);
    else if(decal == 0.)
      t=floor(day*24.*3600.+0.5);
    step=(int)( floor((t-model_start)/info.sampling+0.5)+1 );

    status=clx_archiveread(in1, info, step, buffer1, &time);
    status=clx_archiveread(in2, info, step, buffer2, &time);
    if (init=0) {
      init=1;
      continue;
      }
    printf("t= %6.1f hours frame= %d %s\n",time/3600.0,step,sgetnewdate(info.reference,time));
    for (j=0;j<nndes;j++) {
      H=info.mesh.vertices[j].h;
      h   = buffer3[0][j]-buffer1[0][j];
      p   = buffer4[0][j]-buffer2[0][j];
      ibd=h+p;
      avg[0][j]+=h;       /*elevation*/
      rms[0][j]+=h*h;
      avg[3][j]+=p;       /*air pressure*/
      rms[3][j]+=p*p;
      avg[4][j]+=ibd;     /*inverted barometer departure*/
      rms[4][j]+=ibd*ibd;
      }
    swap=buffer3;
    buffer3=buffer1;
    buffer1=swap;
    swap=buffer4;
    buffer4=buffer2;
    buffer2=swap;
    count++;
    fclose(in1);
    fclose(in2);
    }

    for (i=0;i<nndes;i++) {
      for(k=0;k<nstat;k++) avg[k][i] /=count;
      for(k=0;k<nstat;k++) {
        rms[k][i] /=count;
        rms[k][i] = sqrt(rms[k][i]-avg[k][i]*avg[k][i]);
        }
      };

 /*---------------- save the statistic under bmg format -------------*/
    header.ni=grid.nx;
    header.nj=grid.ny;
    header.xmin=grid.xmin;
    header.ymin=grid.ymin;
    header.dx=grid.dx;
    header.dy=grid.dy;
    header.nk=1;
    header.nt=2;
    header.spec=mask;
    header.levels=(float *)malloc(sizeof(float));
    header.levels[0]=0.;
  
    strcpy(header.comment[0],"Mean and RMS");
    strcpy(header.comment[1],"no comment 1");
    strcpy(header.comment[2],"sea level 10 days variability");
    strcpy(header.comment[3],"no comment 2");

    header.nd=1;
    if(decal !=0.)
      sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-h10d-decal%4.2lf.bmg",zone,start.month,start.year,final.month,final.year,decal);
    else
      sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-h10d.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[0],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[0],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 2, grid, buf,9999.999, mask);

    strcpy(header.comment[0],"Mean and RMS");
    strcpy(header.comment[1],"no comment 1");
    strcpy(header.comment[2],"reduced air pressure 10 days variability");
    strcpy(header.comment[3],"no comment 2");

    header.nd=1;

    if(decal !=0.)
      sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-p10d-decal%4.2lf.bmg",zone,start.month,start.year,final.month,final.year,decal);
    else
      sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-p10d.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[3],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[3],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 2, grid, buf,9999.999, mask);

    strcpy(header.comment[0],"Mean and RMS");
    strcpy(header.comment[1],"no comment 1");
    strcpy(header.comment[2],"IB departure 10 days variability");
    strcpy(header.comment[3],"no comment 2");

    header.nd=1;

    if(decal !=0.)
      sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-ibd10d-decal%4.2lf.bmg",zone,start.month,start.year,final.month,final.year,decal);
    else
      sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-ibd10d.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[4],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[4],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 2, grid, buf,9999.999, mask);


    }

  for(k=0;k<nstat;k++) free(avg[k]);
  for(k=0;k<nstat;k++) free(rms[k]);

  for(k=0;k<3;k++)free(buffer1[k]);
  for(k=0;k<3;k++)free(buffer2[k]);
  for(k=0;k<3;k++)free(buffer3[k]);
  for(k=0;k<3;k++)free(buffer4[k]);

  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
  
}
