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
#include "archive.h"

#define nstat 9

//function de bmg_fe.cpp
extern int fe_filemapr1 (char* input,char* output,grid_t grid,int minstep,int maxstep,meta_archive_t info);
extern int fe_filemapr2 (char* input,char* output,grid_t grid,int minstep,int maxstep,meta_archive_t info);

/*----------------------------------------------------------------------

Compute global statistics for the MOG2D outputs

----------------------------------------------------------------------*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float z,dum[3],dum1[3],dum2[3],velocity[2],elevation,ws[2],pressure;
  float *buffer1,*buffer2;
  double t0,dT,pulsation;
  double *serie[500],a1,p1,a2,p2,zr,zi,d;

  float  h,u,v,p,hu,hv,wsx,wsy,*avg[nstat],*rms[nstat];
  float  count,dummy, hours,ibd;
  date_t start,final;
  double time;
  double *woce_h,*woce_t,woce_mask;
  int woce_n;
  int fmt,i,j,k,l,L,ndum,rstatus,shift,nndes;
  int minstep,maxstep,l1,l2,n,status,nst;
  int *sample,hr,mn,sc;
  int day,month,year,nday,loop,nloop,elapsed;
  size_t size;
  FILE *in6,*out,*sout[500];
  char *hfile,*pfile,*keyword,*zone,*s;
  char *meshfile;
  char file1[256],file2[256],*hroot,*proot;
  char rootname[256]="\0",tmp[256];
  grid_t grid;
  int *list;
  meta_archive_t info;

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
        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

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

        case 'l' :
          s= strdup(argv[n+1]);
          sscanf(s,"%d",&nday);
          n++;
          n++;
          break;

/*
        case 'i' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&day,&month,&year);
          break;

        case 'f' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d/%d/%d",&day,&month,&year);
          break;

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

 /*---------------- read the mesh file (neigbour format) --------------*/

  sprintf(hfile,"%s-%4.4d.%2.2d",hroot,start.year,start.month);
  status= clx_archiveinfo(hfile,&info);
  status=fe_list(&info.mesh);

  nndes=info.mesh.nvtxs;

 /*--------------- Warning, hardcoded time step for output ------------*/

  dT=7200;
  nloop=4392/24/nday;
  
  if (strcmp(zone,"indien")==0) {
  /* Frame indien */
    grid.xmin= +20;
    grid.xmax=+150;
    grid.ymin= -60;
    grid.ymax= +20;
    grid.dx=1.;
    grid.dy=1.;
    }
  else if(strcmp(zone,"medsea")==0) {
    /* Frame medsea */
    grid.xmin=-10;
    grid.xmax=+37.5;
    grid.ymin=+30;
    grid.ymax=+47.5;
    grid.dx=.125;
    grid.dy=.125;
    }
  else if(strcmp(zone,"gibraltar")==0) {
    /* Zoom Gibraltar */
    grid.xmin=-6.5;
    grid.xmax=-4.5;
    grid.ymin=+35.5;
    grid.ymax=+36.5;
    grid.dx=.01;
    grid.dy=.01;
    }
  else if(strcmp(zone,"kerguelen")==0) {
    /* Frame kerguelen */
    grid.xmin=43.;
    grid.xmax=90.;
    grid.ymin=-74.;
    grid.ymax=-35.;
    grid.dx=.1;
    grid.dy=.1;
    }
  
  grid.nx=(int)( 1+(grid.xmax-grid.xmin)/grid.dx );
  grid.ny=(int)( 1+(grid.ymax-grid.ymin)/grid.dy );
  
  grid.modeH=0;
  
  grid.x=(double *) malloc(grid.nx*sizeof(double));
  grid.y=(double *) malloc(grid.ny*sizeof(double));
  
  minstep=0;
  maxstep=0;

  printf("#create bmg maps ...\n");
  for (loop=1;loop<nloop+1;loop++) {

  
    sprintf(file1,"%s%s.%d.%s",rootname,"mean1",loop,"fil");
    l1=strlen(file1);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"meanH",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean elevation gridded file: %s\n",file2);
/*     efc_mapr1(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr1(file1,file2,grid,minstep,maxstep,info);

    sprintf(file2,"%s%s.%d.%s",rootname,"meanV",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean velocity gridded file: %s\n",file2);
/*     efc_mapr2(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr2(file1,file2,grid,minstep,maxstep,info);
  
    sprintf(file1,"%s%s.%d.%s",rootname,"rms1",loop,"fil");
    l1=strlen(file1);

    sprintf(file2,"%s%s.%d.%s",rootname,"rmsH",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean elevation gridded file: %s\n",file2);
/*     efc_mapr1(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr1(file1,file2,grid,minstep,maxstep,info);
 
    sprintf(file2,"%s%s.%d.%s",rootname,"rmsV",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean velocity gridded file: %s\n",file2);
/*     efc_mapr2(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr2(file1,file2,grid,minstep,maxstep,info);
  
    sprintf(file1,"%s%s.%d.%s",rootname,"mean2",loop,"fil");
    l1=strlen(file1);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"meanD",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean depth gridded file: %s\n",file2);
/*     efc_mapr1(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr1(file1,file2,grid,minstep,maxstep,info);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"meanT",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean flux gridded file: %s\n",file2);
/*     efc_mapr2(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr2(file1,file2,grid,minstep,maxstep,info);
  
    sprintf(file1,"%s%s.%d.%s",rootname,"rms2",loop,"fil");
    l1=strlen(file1);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"rmsD",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean depth gridded file: %s\n",file2);
    /*    efc_mapr1(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus);*/
    status=fe_filemapr1(file1,file2,grid,minstep,maxstep,info);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"rmsT",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean flux gridded file: %s\n",file2);
/*     efc_mapr2(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr2(file1,file2,grid,minstep,maxstep,info);

    sprintf(file1,"%s%s.%d.%s",rootname,"mean3",loop,"fil");
    l1=strlen(file1);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"meanP",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean pressure gridded file: %s\n",file2);
/*     efc_mapr1(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr1(file1,file2,grid,minstep,maxstep,info);
   
    sprintf(file2,"%s%s.%d.%s",rootname,"meanW",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean flux gridded file: %s\n",file2);
/*     efc_mapr2(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr2(file1,file2,grid,minstep,maxstep,info);

    sprintf(file1,"%s%s.%d.%s",rootname,"rms3",loop,"fil");
    l1=strlen(file1);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"rmsP",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean pressure gridded file: %s\n",file2);
/*     efc_mapr1(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr1(file1,file2,grid,minstep,maxstep,info);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"rmsW",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean flux gridded file: %s\n",file2);
/*     efc_mapr2(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr2(file1,file2,grid,minstep,maxstep,info);

    sprintf(file1,"%s%s.%d.%s",rootname,"mean4",loop,"fil");
    l1=strlen(file1);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"meanIBD",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean pressure gridded file: %s\n",file2);
/*     efc_mapr1(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr1(file1,file2,grid,minstep,maxstep,info);
  
    sprintf(file1,"%s%s.%d.%s",rootname,"rms4",loop,"fil");
    l1=strlen(file1);
  
    sprintf(file2,"%s%s.%d.%s",rootname,"rmsIBD",loop,"bmg");
    l2=strlen(file2);
    fprintf(stderr,"create a mean pressure gridded file: %s\n",file2);
/*     efc_mapr1(file1,file2,grid,Nnodes,minstep,maxstep,&rstatus); */
    status=fe_filemapr1(file1,file2,grid,minstep,maxstep,info);
  

    }
  fclose(in6);

  free(list);

   __ERR_BASE_LINE__("exiting\n");exit(0);
  
}

