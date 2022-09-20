#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "poc-time.h"
#include "fe.h"
#include "archive.h"
#include "functions.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_intpl_xyt(char *datafile[2], meta_archive_t info, double t, double p, double time,float serie[3], memory_t *memory)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  time : julian days, CNES time origin

-----------------------------------------------------------------------------*/
{
  float  z,z1,z2;
  double beta[3];
  float  count,start,end;
  double time1,time2,sampling=0.0,t0;
  int nodes[3];
  int i,j,k,node;
  int status;
  int nndes,elt,read;
  date_t reference;


  FILE *in1=NULL, *in2=NULL;

  reference=info.reference;

  t0=(julian_day(reference.month,reference.day,reference.year)-CNES0jd)*24.*3600.+reference.second;

  start=t0+info.start;

  status=fe_beta(info.mesh, t, p,&elt,nodes,beta);
  printf("%lf %lf \n", t,p);
  printf("%d %d %d %d %lf %lf %lf \n", elt, nodes[0],nodes[1],nodes[2],beta[0],beta[1],beta[2]);
  if(status!=0) return(-1);

  nndes=info.mesh.nvtxs;

  read=1;

  if(memory->file==NULL) goto skip;

  if(strcmp(memory->file,datafile[0])==0) {
    i=(int) floor ((t-start)/info.sampling);
    if(i==memory->frame) read=0;
    else read=1;
    }

skip:
   if (read) {
     in1=fopen(datafile[0],"r");
     if(memory->file!=NULL) free(memory->file);
     memory->file=strdup(datafile[0]);
     }

    t=time*24.*3600.;
    i=(int) floor ((t-start)/info.sampling);
    if (read) {
      status=clx_archiveread(in1,info,i+1,memory->h1,&time1);
      if(status!=0) {
        goto error;
        }
      if(i+2<=info.nframe) {
        status=clx_archiveread(in1,info,i+2,memory->h2,&time2);
        }
      else {
        in2=fopen(datafile[1],"r");
        status=clx_archiveread(in2,info,1,memory->h2,&time2);
        fclose(in2);
        }
      if(status!=0) {
        goto error;
        }
      memory->frame=i;
      }
    printf("t= %6.1f hours status= %d %s\n",time1/3600.0,status, sgetnewdate(info.reference,time1));
    printf("t= %6.1f hours status= %d %s\n",time2/3600.0,status, sgetnewdate(info.reference,time2));
    t=t-t0;
    for (j=0;j<3;j++) {
      z1=0;
      z2=0;
      for (i=0;i<3;i++) {
//         z1+=memory->h1[j][nodes[i]-1]*beta[i];
//         z2+=memory->h2[j][nodes[i]-1]*beta[i];
        z1+=memory->h1[j][nodes[i]]*beta[i];
        z2+=memory->h2[j][nodes[i]]*beta[i];
        }
      serie[j] = ((time2-t)*z1+(t-time1)*z2)/info.sampling;
      }

  if (in1!=0) fclose(in1);
  return(0);

  error:

  if (in1!=0) fclose(in1);
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dlon,dlat;
  float serie[3];
  char point_name[50];
  char a;
  int Mes,countnbdata;
  char line[1000];

  double lon,lat;
  float  h,p,z;
  float  dummy, hours,start,end;
  double time,t;
  double time1,time2,sampling=0.0;
  int origine;
  int i,j,k,l,L,count;
  int n,status,nndes;
  int status1,status2,status3,status4;
  int year,month,day,hour;
  date_t reference,actual,current;

  FILE *in,*out;
  char *rootname_h=NULL,*rootname_p=NULL;
  char *option,*keyword,*s,*datafile_h[2],*datafile_p[2],*input,*output;
  char *default_dir=0,*pathname=0;
  meta_archive_t info;
  memory_t memory_h,memory_p;
  char *check;

  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'h' :
          rootname_h= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          rootname_p= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*---------------------------------------------------------------------
        pathname for default output directory */
//         case 'd' :
//           default_dir= strdup(argv[n+1]);
//           n++;
//           n++;
//           break;

/*---------------------------------------------------------------------
        pathname for default archive directory */
        case 'r' :
          pathname= strdup(argv[n+1]);
          n++;
          n++;
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


  exitIfNull(
    in=fopen(input,"r")
    );

//   if(default_dir==NULL) {
//     default_dir=strdup(".");
//     printf("use <.> as root name for default output directorys\n");
//     }

  exitIfNull(
    out=fopen(output,"w")
    );

  if(pathname==NULL) {
    pathname=strdup(".");
    printf("use <.> as root name for default archives directorys\n");
    }

  info.mesh.triangles=NULL;

  memory_h.activated=0;
  memory_p.activated=0;
  memory_h.file=NULL;
  memory_p.file=NULL;

  // Lecture du nombre de points de mesure

  do {a=fgetc(in);} while ( a!= ':');
  fscanf(in, "%d",&Mes);
  do {a=fgetc(in);} while ( a!= '\n');

  printf("Nbre de mesures : %d\n",Mes);
  
  current.month=-1;

  // Lecture des données

  for (i=0; i<Mes;i++) {
      check=fgets(line,sizeof(line),in);
      if(check==0) goto finished;
      next :
      sscanf(line, " %lf %lf %lf %s", &lon,&lat,&time,&point_name);

      actual=poctime_getdatecnes(time*24,'h');
/*     printf("%s %f %f\n", sgetcnesdate(time*24),lon,lat); */

      datafile_h[0]=(char *) malloc(strlen(pathname)+strlen(rootname_h)+64);
      datafile_h[1]=(char *) malloc(strlen(pathname)+strlen(rootname_h)+64);

      datafile_p[0]=(char *) malloc(strlen(pathname)+strlen(rootname_p)+64);
      datafile_p[1]=(char *) malloc(strlen(pathname)+strlen(rootname_p)+64);

      sprintf(datafile_h[0],"%s/%s-%4.4d.%2.2d",pathname,rootname_h,actual.year,actual.month);
      sprintf(datafile_p[0],"%s/%s-%4.4d.%2.2d",pathname,rootname_p,actual.year,actual.month);

      if(actual.month!=current.month) {
        printf("opening %s %s\n",datafile_h[0],datafile_p[0]);
        status1=clx_archiveinfo(datafile_h[0],&info);
        if(status1==0) {
          nndes=info.mesh.nvtxs;
          if(memory_h.activated==0) {
            memory_h.h1=smatrix(0,2,0,nndes-1);
            memory_h.h2=smatrix(0,2,0,nndes-1);
            memory_h.activated=1;
            }
          if(memory_p.activated==0) {
            memory_p.h1=smatrix(0,2,0,nndes-1);
            memory_p.h2=smatrix(0,2,0,nndes-1);
            memory_p.activated=1;
            }
          reference=info.reference;
          status2=fe_list(&info.mesh);
          if(status2!=0) goto error;
          current=actual;
          }
        else goto next;
        }
      date_t next=actual=poctime_getdatecnes(time*24+info.sampling/3600.,'h');
      sprintf(datafile_h[1],"%s/%s-%4.4d.%2.2d",pathname,rootname_h,next.year,next.month);
      sprintf(datafile_p[1],"%s/%s-%4.4d.%2.2d",pathname,rootname_p,next.year,next.month);
      
      status3=fe_intpl_xyt(datafile_h,info,lon,lat,time,serie,&memory_h);
      if(status3==0) {
        h=serie[0];
        status4=fe_intpl_xyt(datafile_p,info,lon,lat,time,serie,&memory_p);
        if(status4==0) {
          p=serie[0];
          fprintf(out,"%lf %lf %lf %s %f %f \n",lon,lat,time,point_name,h,p);
/*        fprintf(out,"%lf %lf %lf %f %f %f \n",t,dlat,dlon,z,h,p); */
          }
        }
    }

finished:
  fclose(in);
  fclose(out);
  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  fprintf(stderr,"%s -computation aborted ^^^^^^^^^\n",argv[0]);
  fclose(in);
  fclose(out);
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}

