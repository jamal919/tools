#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
 
#include "tools-structures.h"
    
#include "rutin.h"     /*  rutin.h contains common utility routines  */
#include "archive.h"
#include "fe.h"
#include "poc-time.h"

date_t reference;
double t0; /*attention au traitement de simulation n'ayant pas la meme reference*/
list_t list;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_intpl_xyt1_old(char * datafile,char * datafile2, meta_archive_t info, double t, double p, double time,float serie[3], memory_t *memory,float start) __attribute__((warning("This calls deprecated fe_beta_inlist")));
int fe_intpl_xyt1_old(char * datafile,char * datafile2, meta_archive_t info, double t, double p, double time,float serie[3], memory_t *memory,float start)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  beta[3];
  float  z,z1,z2;
  float  count,end;
  double time1,time2,sampling=0.0,tt;
  int origine, nodes[3];
  int i,j,k,l,node;
  int kk,ll;
  int status;
  int nndes,elt,read;
  meta_archive_t info2;
  FILE *in,*in2;
    
  k=(int)( floor(t/10.) );
  l=(int)( floor(p/10.)+9 );
  if(k<0)   k+=36;
  if(k>=36) k-=36;
  
  status=fe_beta_inlist(info.mesh,list.elements[k][l],t, p,&elt,nodes,beta);
  if(status!=0) {
/*
    if((k==33)&&(l>2)) {
      printf("%d %d %f %f %lf %d\n",k,l,t,p,time,status);
      for (kk=k-1;kk<=k+1;kk++) for (ll=l-1;ll<=l+1;ll++) {
  status=fe_beta_inlist(info.mesh,list[kk][ll],nlisted[kk][ll],t, p,&elt,nodes,beta);
      printf("%d %d %d\n",kk,ll,status);
        }
    }
*/
    return(-1);
    }

  nndes=info.mesh.nvtxs;

  read=1;

  if(memory->file==NULL) {
    goto skip;
    }

  if(strcmp(memory->file,datafile)==0) {
    tt=time*24.*3600.;
    i=(int)( (tt-start)/info.sampling );
    if(i==memory->frame) read=0;
    else read=1;
    }
   
    
skip:
   if (read) {
     in=fopen(datafile,"r");
     if(in == NULL) printf("error while opening the datafile %s\n",datafile);
     if(memory->file!=NULL) free(memory->file);
     memory->file=strdup(datafile);
     }
     
  tt=time*24.*3600.;
  i=(int)( (tt-start)/info.sampling );
  if((i>=info.nframe-1) && (read)) {
    i=info.nframe;
    status=clx_archiveread(in,info,i,memory->h1,&memory->time1);
    if(status!=0) goto error;
    memory->frame=i-1;
    
    status=clx_archiveinfo(datafile2,&info2);
    in2=fopen(datafile2,"r");
    if(in2 == NULL) printf("error while opening the datafile2 %s\n",datafile2);
    i=1;
    status=clx_archiveread(in2,info2,i,memory->h2,&memory->time2);
    if(status!=0) goto error;
    archive_freeinfo(&info2);
    fclose(in2);
    }
  else if((i<info.nframe-1) &&  (read)) {
    status=clx_archiveread(in,info,i+1,memory->h1,&memory->time1);
    if(status!=0) goto error;
    status=clx_archiveread(in,info,i+2,memory->h2,&memory->time2);
    if(status!=0) goto error;
    memory->frame=i;
/*    printf("#1 read for t= %6.1f hours status= %d %s topex= %lf\n",memory->time1/3600.0,status, sgetnewdate(info.reference,memory->time1),(time*24-2922*24)*3600);
    printf("#2 read for t= %6.1f hours status= %d %s\n",memory->time2/3600.0,status, sgetnewdate(info.reference,memory->time2));
*/
    }
/* printf("t= %6.1f hours status= %d %s\n",memory->time1/3600.0,status, sgetnewdate(info.reference,memory->time1));
   printf("t= %6.1f hours status= %d %s\n",memory->time2/3600.0,status, sgetnewdate(info.reference,memory->time2)); */
  tt=tt-t0;
  for (j=0;j<3;j++) {
    z1=0;
    z2=0;
    for (i=0;i<3;i++) {
      z1+=memory->h1[j][nodes[i]-1]*beta[i];
      z2+=memory->h2[j][nodes[i]-1]*beta[i];
      }
    serie[j] = ((memory->time2-tt)*z1+(tt-memory->time1)*z2)/info.sampling;
    }
 
  if (read) fclose(in);
  return(0);

  error:

  fclose(in);
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_intpl_xyt1(FILE *in, char *datafile, meta_archive_t info, double t, double p, double time,float serie[3], memory_t *memory);
int fe_intpl_xyt1(FILE *in, char *datafile, meta_archive_t info, double t, double p, double time,float serie[3], memory_t *memory)
//not used
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  z,z1,z2;
  double beta[3];
  float  count,start,end;
  double time1,time2,sampling=0.0,t0,tt;
  int origine, nodes[3];
  int i,j,k,node;
  int status;
  int nndes,elt,read;
  date_t reference;

  reference=info.reference;

  origine=julian_day(1,1,1950);
  t0=(julian_day(reference.month,reference.day,reference.year)-origine)*24.*3600.+reference.second;

  start=t0+info.start;

  status=fe_beta(info.mesh, t, p,&elt,nodes,beta);
  if(status!=0) return(-1);

  nndes=info.mesh.nvtxs;
  read=1;
   
  if (in == NULL) {
     in=fopen(datafile,"r");
     if(in == NULL) printf("error while opening the datafile %s\n",datafile);
     if(memory->file!=NULL) free(memory->file);
     memory->file=strdup(datafile);
     }
     
  tt=time*24.*3600.;
  i=(int)( (tt-start)/info.sampling );
  if (read) {
    status=clx_archiveread(in,info,i+1,memory->h1,&memory->time1);
    if(status!=0) goto error;
    status=clx_archiveread(in,info,i+2,memory->h2,&memory->time2);
    if(status!=0) goto error;
    memory->frame=i;
    }
/*   printf("t= %6.1f hours status= %d %s\n",memory->time1/3600.0,status, sgetnewdate(info.reference,memory->time1));
   printf("t= %6.1f hours status= %d %s\n",memory->time2/3600.0,status, sgetnewdate(info.reference,memory->time2)); */
  tt=tt-t0;
  for (j=0;j<3;j++) {
    z1=0;
    z2=0;
    for (i=0;i<3;i++) {
      z1+=memory->h1[j][nodes[i]-1]*beta[i];
      z2+=memory->h2[j][nodes[i]-1]*beta[i];
      }
    serie[j] = ((memory->time2-tt)*z1+(tt-memory->time1)*z2)/info.sampling;
    }
 
  return(0);

  error:
   return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double dlon,dlat;
  float serie[3];
  
  double  lon,lat;
  float  z;
  float  dummy, hours,starter,end;
  double time,dum,second;
  double time1,time2,sampling=0.0;
  int origine,cycle,nitems,h,p;
  int i,j,k,l,L,count,minute,flag;
  int n,status,nndes,t,mask=999999;
  int year,month,day,hour;
  date_t reference,actual,current,current2,start_date;

  FILE *in1,*in2,*in3,*out,*inh=NULL,*inp=NULL;
  char *rootname_h=NULL,*rootname_p=NULL;
  char *option,*keyword,*s,*datafile_h=NULL,*datafile_p=NULL,*input1,*input2,*output;
  char *datafile2_h=NULL,*datafile2_p=NULL;
  meta_archive_t info;
  memory_t memory_h,memory_p;
  
  extern int fe_nlistmax;
 
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
        case 'c' :
          s= strdup(argv[n+1]);
          n++;
          n++;
          sscanf(s,"%d",&cycle);
          break;

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

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        break;
      }
      free(keyword);
    }

  input1=(char *)malloc(1024);
  sprintf(input1,"cyclin_topex_%3.3d",cycle);

  output=(char *)malloc(1024);
  sprintf(output,"MOG2D_topex_%3.3d",cycle);

  in1=fopen(input1,"r");

  out=fopen(output,"w");

  info.mesh.triangles=NULL;

  memory_h.activated=0;
  memory_p.activated=0;
  memory_h.file=NULL;
  memory_p.file=NULL;
  
  list.init(36,18);
    
  do{
    next:
    /* lat et lon en degres; t en secondes depuis 1/1/1958 */
    /* nb jours de 1950 -> 1958 = 2922                     */
    nitems=fscanf(in1,"%f %f %ld",&lat,&lon,&t);
    /*printf("%f %f %ld\n",lat,lon,t);*/
    if(nitems!=3) goto end;

    time=t/(3600.)+2922.*24.;

    actual=poctime_getdatecnes(time,'h');
    /*printf("#tps read :%s , JJ1950 = %lf, sec T/P = %ld\n", sgetcnesdate(time),time/24.,t);*/
    
    if(datafile_h==NULL) datafile_h=(char *)malloc(strlen(rootname_h)+64);
    if(datafile_p==NULL) datafile_p=(char *)malloc(strlen(rootname_p)+64);
    
    sprintf(datafile_h,"%s-%4.4d.%2.2d.reduced-6",rootname_h,actual.year,actual.month);
    sprintf(datafile_p,"%s-%4.4d.%2.2d.reduced-6",rootname_p,actual.year,actual.month);

    if(actual.month!=current.month) {
      printf("opening %s %s\n",datafile_h,datafile_p);
      archive_freeinfo(&info);
      status=clx_archiveinfo(datafile_h,&info);
      if(status!=0) goto next;
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
      origine=julian_day(1,1,1950);
      t0=(julian_day(reference.month,reference.day,reference.year)-origine)*24.*3600.+reference.second;
      starter=t0+info.start;
      status=fe_list(&info.mesh);
      if(status!=0) goto error;
      grid_t RIEN; //modif Thierry
      fe_CreateTriangleList(RIEN,info.mesh,&list);
      current=actual;
      }
    else
      {
      if(actual.second >= 64800.) {
        start_date=poctime_getdatecnes(time+6.,'h');
        if((start_date.month != actual.month) && (start_date.month != current2.month)) {
          if(datafile2_h==NULL) datafile2_h=(char *)malloc(strlen(rootname_h)+64);
          if(datafile2_p==NULL) datafile2_p=(char *)malloc(strlen(rootname_p)+64);
          sprintf(datafile2_h,"%s-%4.4d.%2.2d.reduced-6",rootname_h,start_date.year,start_date.month);
          sprintf(datafile2_p,"%s-%4.4d.%2.2d.reduced-6",rootname_p,start_date.year,start_date.month);
          current2=start_date;
          }
        }
      }

    time=time/24.;
    status=fe_intpl_xyt1_old(datafile_h,datafile2_h,info,lon,lat,time,serie,&memory_h,starter);
/*  hauteur modele */
    if(status!=0) h=mask;
    else h=(int)(serie[0]*1000);

    status=fe_intpl_xyt1_old(datafile_p,datafile2_p,info,lon,lat,time,serie,&memory_p,starter);
/*  hauteur IB */
    if(status!=0) p=mask;
    else p=(int)(-1.*serie[0]*1000);

    fprintf(out,"%6d %6d\n",h,p);
/*    printf("%6d %6d\n",h,p);*/
    
    } while (!feof(in1));

 end:
  fclose(in1);
  fclose(out);
  __ERR_BASE_LINE__("exiting\n");exit(0);

 error:
  fprintf(stderr,"%s -computation aborted ^^^^^^^^^\n",argv[0]);
  fclose(in1);
  fclose(out);
  __ERR_BASE_LINE__("exiting\n");exit(-1);
}

