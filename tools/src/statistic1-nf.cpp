
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

#include "tools-define.h"
#include "tools-structures.h"

#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "map.h"
#include "fe.h"
#include "poc-time.h"
#include "archive.h"
#include "bmg.h"
#include "functions.h"


#define nstat 10

/*----------------------------------------------------------------------

Compute global statistics for the MOG2D outputs

----------------------------------------------------------------------*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_zonegrid_fatal(char *zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;
  int zone_initialised=0;

  /* BISCAY grid */
  if(strcmp(zone,"BISCAY")==0) {
  map_set2Dgrid(&grid,-13.0,+42.0,+0.0,+53.0,0.1,0.1);
  grid.nx  = 231;
  grid.ny  = 121;
  grid.modeH=0;
  zone_initialised=1;
  }
  
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
  grid.xmin= +43.0 ;
  grid.ymin= -73.0;
  grid.dx  =   0.1;
  grid.dy  =   0.1;
  grid.nx  = 471;
  grid.ny  = 381;
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
      map_set2Dgrid(&grid,0.,-90.,360.,90.,1,1);
      grid.nx=361;
      grid.ny=181;
      grid.modeH=0;
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
  grid.modeH=0;/* ??????? LR, add 11/07/05 */
  zone_initialised=1;
  }
  
/* mesdsea grid, normal resolution */
  if(strcmp(zone,"medsea")==0) {
  map_set2Dgrid(&grid,-10.0,27.5,40.0,47.5,0.25,0.25);
  grid.nx  = 201;
  grid.ny  = 81;
  grid.modeH=0;/* ??????? LR, add 18/05/04 */
  zone_initialised=1;
  }

/* caspian grid, normal resolution, LR, add 23/06/05 */
 if(strcmp(zone,"caspian")==0) {
    /* Frame global */
    map_set2Dgrid(&grid,46.0,36.0,55.0,47.5,0.05,0.05);
    grid.modeH=0;
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
  float *buffer1[3],*buffer2[3],*buf,*bufx,*bufy;

  float  H,h,u,v,p,hu,hv,wsx,wsy,hwsx,hwsy,*avg[nstat],*rms[nstat];
  float  count, hours,mask=1.e+10,ibd;
  double time,t;
  int step,nndes;
  int fmt,i,j,k,l,L;
  int l1,l2,n,status;
  int day,month,year,nday,loop,elapsed;
  int *element;
  bmgheader_t header;

  size_t size;
  FILE *in1,*in2,*out,*fout;
  char hfile[1024],pfile[1024],*hroot,*proot,*keyword,*zone,*s;

  char file1[256],file2[256],mean_P[1024],mean_SL[1024];
  char rootname[256]="\0",tmp[256],output[1024],input[1024];
  grid_t grid;
  meta_archive_t info;

  date_t start,final;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
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

  if(hroot==NULL) {
    hroot= strdup("analysis");
    printf("use <analysis> as root name for sea state archive files\n");
    }
  if(proot==NULL) {
    proot= strdup("forcing");
    printf("use <forcing> as root name for forcing archive files\n");
    }
  printf("rootname: %s \n",rootname);

 /*---------------- read the mesh file (neigbour format) --------------*/

  sprintf(hfile,"%s-%4.4d.%2.2d",hroot,start.year,start.month);
  status= clx_archiveinfo(hfile,&info);
  status=fe_list(&info.mesh);

  nndes=info.mesh.nvtxs;

 /*-------------------- allocate memory for vectors -------------------*/

  for(k=0;k<nstat;k++) avg[k]=(float *)malloc(nndes*sizeof(float));
  for(k=0;k<nstat;k++) rms[k]=(float *)malloc(nndes*sizeof(float));

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
    }

  for(loop=0;loop<1;loop++) {
    for (i=0;i<nndes;i++) {
      for(k=0;k<nstat;k++) avg[k][i]=0.;
      for(k=0;k<nstat;k++) rms[k][i]=0.;
      };

  grid=get_zonegrid_fatal(zone);
  status=map_completegridaxis_2(&grid);

  element=fe_scan_elements(info.mesh,grid,0);
  if (element == NULL) {__ERR_BASE_LINE__("exiting\n");exit(8);}

  buf=(float *)malloc(grid.nx*grid.ny*sizeof(float));
  if (buf == NULL) {__ERR_BASE_LINE__("exiting\n");exit(9);}

  count=0;
  for (year=start.year;year<=final.year;year++) {
    for (month=1;month<=12;month++) {
      if ((year==start.year) && (month<start.month)) continue;
      if ((year==final.year) && (month>final.month)) break ;

      sprintf(hfile,"%s-%4.4d.%2.2d",hroot,year,month);
      sprintf(pfile,"%s-%4.4d.%2.2d",proot,year,month);

      in1=fopen(hfile,"r");
      in2=fopen(pfile,"r");

      printf("read %s & %s ...\n",hfile,pfile);

 /*---------------- read the output files (MOG2D format) --------------*/

      status= clx_archiveinfo(hfile,&info);

      for(step=1;step<=info.nframe;step++) {
        status=clx_archiveread(in1, info, step, buffer1, &time);
        status=clx_archiveread(in2, info, step, buffer2, &time);
        if(step%25 ==0) printf("t= %6.1f hours frame= %d %s\n",time/3600.0,step,sgetnewdate(info.reference,time));
        for (j=0;j<nndes;j++) {
          H=info.mesh.vertices[j].h;
          h   = buffer1[0][j];
          p   = buffer2[0][j];
          u   = buffer1[1][j];
          wsx = buffer2[1][j];
          v   = buffer1[2][j];
          wsy = buffer2[2][j];
          ibd=h+p;
          hu=(H+h)*u;
          hv=(H+h)*v;
          hwsx=h*wsx;
          hwsy=h*wsy;
          avg[0][j]+=h;       /*elevation*/
          rms[0][j]+=h*h;
          avg[1][j]+=u;       /*zonal current*/
          rms[1][j]+=u*u;
          avg[2][j]+=v;       /*meridional current*/
          rms[2][j]+=v*v;
          avg[3][j]+=p;       /*air pressure*/
          rms[3][j]+=p*p;
          avg[4][j]+=ibd;     /*inverted barometer departure*/
          rms[4][j]+=ibd*ibd;
          avg[5][j]+=hu;      /*zonal transport*/
          rms[5][j]+=hu*hu;
          avg[6][j]+=hv;      /*meridional transport*/
          rms[6][j]+=hv*hv;
          avg[7][j]+=hwsx;    /*zonal depth-integrated wind stress*/
          rms[7][j]+=hwsx*hwsx;
          avg[8][j]+=hwsy;    /*meridional depth-integrated wind stress*/
          rms[8][j]+=hwsy*hwsy;
          avg[9][j]+=u*u+v*v;       /*mean kinetic energy*/
          rms[9][j]+=(u*u+v*v)*(u*u+v*v);
          }
        }
      count+=info.nframe;
      fclose(in1);
      fclose(in2);
      }
    }

 /*---------------- save the statistic under MOG2D format -------------*/

    for (i=0;i<nndes;i++) {
      for(k=0;k<nstat;k++) avg[k][i] /=count;
      for(k=0;k<nstat;k++) {
        rms[k][i] /=count;
        rms[k][i] = sqrt(rms[k][i]-avg[k][i]*avg[k][i]);
        }
      };

 /*---------------- save the mean pressure under ascii format -------------*/
 /*---------------- rajout�le 1/10/2002 par LC               -------------*/
 /**/
    sprintf(mean_P,"MeanPressure_%2.2d_%4d-%2.2d_%4d.s2r",start.month,start.year,final.month,final.year);
    fout=fopen(mean_P,"w");
    fprintf(fout,"Mean Atmospheric Pressure for the period 92-2002 \n");
    fprintf(fout,"in mbars \n");
    for (i=0;i<nndes;i++) {
      fprintf(fout,"%d \t %f \n",i+1,avg[3][i]);
      }
    fclose(fout);

    sprintf(mean_SL,"MeanSlevel_%2.2d_%4d-%2.2d_%4d.s2r",start.month,start.year,final.month,final.year);
    fout=fopen(mean_SL,"w");
    fprintf(fout,"Mean Sea Level for the period 92-2002 \n");
    fprintf(fout,"in cm \n");
    for (i=0;i<nndes;i++) {
      fprintf(fout,"%d \t %f \n",i+1,avg[0][i]);
      }
    fclose(fout);
    goto end;
/**/
 
 /*---------------- save the statistic under MOG2D format -------------*/
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
    strcpy(header.comment[2],"sea level");
    strcpy(header.comment[3],"no comment 2");

    header.nd=1;

    sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-h.bmg",zone,start.month,start.year,final.month,final.year);
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
    strcpy(header.comment[2],"reduced air pressure");
    strcpy(header.comment[3],"no comment 2");

    header.nd=1;

    sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-p.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[3],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[3],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 2, grid, buf, 9999.999,mask);

    strcpy(header.comment[0],"Mean and RMS");
    strcpy(header.comment[1],"no comment 1");
    strcpy(header.comment[2],"IB departure");
    strcpy(header.comment[3],"no comment 2");

    header.nd=1;

    sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-ibd.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[4],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[4],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 2, grid, buf,9999.999, mask);

    strcpy(header.comment[0],"Mean and RMS");
    strcpy(header.comment[1],"no comment 1");
    strcpy(header.comment[2],"Barotropic currents");
    strcpy(header.comment[3],"no comment 2");

    header.nd=2;

    sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-u.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[1],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,avg[2],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 2, 1, grid, buf, 9999.999,mask);
    status=fe_map(info.mesh,rms[1],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 2, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[2],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 2, 2, grid, buf,9999.999, mask);

    strcpy(header.comment[0],"Mean and RMS");
    strcpy(header.comment[1],"no comment 1");
    strcpy(header.comment[2],"Barotropic transport");
    strcpy(header.comment[3],"no comment 2");

    header.nd=2;

    sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-t.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[5],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,avg[6],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 2, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[5],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 2, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[6],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 2, 2, grid, buf,9999.999, mask);

    strcpy(header.comment[0],"Mean and RMS");
    strcpy(header.comment[1],"no comment 1");
    strcpy(header.comment[2],"Depth integrated wind stress");
    strcpy(header.comment[3],"no comment 2");

    header.nd=2;

    sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-ws.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[7],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,avg[8],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 2, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[7],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 2, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[8],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 2, 2, grid, buf, 9999.999,mask);

    strcpy(header.comment[0],"Mean and RMS");
    strcpy(header.comment[1],"no comment 1");
    strcpy(header.comment[2],"Kinetic energy");
    strcpy(header.comment[3],"no comment 2");

    header.nd=1;

    sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-ke.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[9],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf,9999.999, mask);
    status=fe_map(info.mesh,rms[9],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 2, grid, buf,9999.999, mask);

    for (i=0;i<nndes;i++) {
      avg[9][i] = sqrt(avg[9][i]);
      };

    strcpy(header.comment[0],"Mean");
    strcpy(header.comment[1],"no comment 1");
    strcpy(header.comment[2],"current amplitude");
    strcpy(header.comment[3],"no comment 2");

    header.nd=1;
    header.nt=1;

    sprintf(output,"%s-%2.2d.%4d-%2.2d.%4d-u0.bmg",zone,start.month,start.year,final.month,final.year);
    out=fopen(output,"wb");
    status=bmg_writeheader (out, header);
    if(status!=0) goto error;
    fclose(out);

    status=fe_map(info.mesh,avg[9],grid,element,buf,mask);
    status=bmg_saver1 (output, 1, 1, 1, grid, buf, 9999.999,mask);

    }

  end:
  fclose(in1);
  fclose(in2);

  for(k=0;k<nstat;k++) free(avg[k]);
  for(k=0;k<nstat;k++) free(rms[k]);

  for(k=0;k<3;k++)free(buffer1[k]);
  for(k=0;k<3;k++)free(buffer2[k]);

  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  __ERR_BASE_LINE__("exiting\n");exit(-1);
  
}
