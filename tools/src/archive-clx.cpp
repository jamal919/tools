#include <config.h>
   
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "constants.h"
 
#include "fe.h"
#include "archive.def"
#include "archive.h"
#include "poc-time.h"
#include "parallel.h"


//define de swap function to be homegenous with the Sun Sparc binary format
//-------------------------------------------------------------------------
#if NEED_SWAP == 1
#include "swap.h"
#endif
//-------------------------------------------------------------------------

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int identify_BINARY(const char *file, const char *name, int *id, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
  if(strcmp(name,"elevation")==0) {
    *id=0;
    }
  if(strcmp(name,"ubar")==0) {
    *id=1;
    }
  if(strcmp(name,"vbar")==0) {
    *id=2;
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int Read_TugoBinary(const char *datafile, double time, float **buffer, double *actual_time, double *next)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int step, status, d1, d2;
  double start,shift;

  FILE *in;
  date_t actual;
  meta_archive_t info;
//   bool spokesman=(gCPU_ID==gCPU_MASTER);
  
//   getmodeldate(time, &actual);
//   char *sdate=sgetmodeldate(time);
//   if(spokesman) printf("#obc_ReadBinary at %s in %s\n", sgetmodeldate(time), datafile);
//   free(sdate);
  
  status = archive_info(datafile, TUGO_UG_BINARY_FORMAT, &info);
  if(status != 0) {
    printf("error while reading header in datafile: %s\n",datafile);
    goto error;
    }

/**---------------------------------------------------------------------
  model time reference */
  d1 = julian_day(t_reference.month, t_reference.day, t_reference.year);
  
/**---------------------------------------------------------------------
  external archive time reference */
  d2 = julian_day(info.reference.month, info.reference.day,info.reference.year);

/**---------------------------------------------------------------------
  shift between model and archive time: archive time=model time -start */
  shift = (d2 - d1) * 24 * 3600 - t_reference.second;
  start = info.start + shift;

  step = 1 + (int) ((time - start) / info.sampling);
  in = fopen(datafile, "r");
  if(in == NULL) {
    printf("error while opening the datafile %s\n", datafile);
    goto error;
    }

  status = clx_archiveread(in, info, step, buffer, actual_time);
  if(status != 0) {
    printf("error while reading the datafile at step =%d %s\n", step,datafile);
//     if(step < 1)
//       printf("model start before archive content : %s\n", sgetmodeldate(start));
    goto error;
    }

  fclose(in);

/**---------------------------------------------------------------------
  convert archive time in model time */
  *actual_time += shift;
  *next = *actual_time +info.sampling;
     
  info.destroy();

  return (status);

error:
  status = -1;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//   int stream_UGxBINARY(UGfield_t<double> & field, string units, double t)
  int stream_UGxBINARY(UGfield_t<double> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   using namespace Keywords ;
  const string KEY_MKS = "MKS" ;
  const string KEY_CGS = "CGS" ;
  
  string units = KEY_MKS;
  
  int status, k;
  double time,next;
  float *buffer[3];
  filestream_t<float> *stream;
  date_t actual, reference;
  double scale;
  stream=(filestream_t<float> *) field.stream;
  if(stream==0) return(-1);
  
  stream->check(t);
  
//   if(external_obc == 0)
//     return (0);

/*------------------------------------------------------------------------------
  scale factor for backward compatibility*/
  if (units == KEY_MKS) {
    scale=1.0;
    }
  if (units == KEY_CGS) {
    scale=1.0e-02;
    }

  for(k = 0; k < 3; k++)
    buffer[k] = new float[field.descriptor->nnodes];
  
  status = Read_TugoBinary(stream->filename.c_str(), t, buffer, &time, &next);
  if(status!=0) return (-1);
  
  field.time=time;
  field.next=next;
  for(int n=0;n<field.descriptor->nnodes;n++) field.x[n]=scale*buffer[stream->id][n];

  field.units=poc_strdup("UNDOCUMENTED");

  for(k = 0; k < 3; k++)
    deletep(&buffer[k]);
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void getnewdate(date_t reference, double t,date_t *actual)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  long N,nday,N1,N2;
  float second;

/* *----------------------------------------------------------------------------
  Get the calendary date corresponding to absolute time t_reference+t */

/*   printf("t_reference is: %f %2d/%2d/%4d \n", */
/*     t_reference.second,t_reference.day,t_reference.month,t_reference.year); */
  N1=julian_day(reference);
  N2=julian_day(1,1,1950);

  N=N1-N2;

/*   printf("N is: %d %d %d \n",N,N1,N2); */
  nday=N+ int( (reference.second+t)/(24*3600.) );

/*   printf("nday is: %d \n",nday); */
  calendary(nday,actual);
  second=fmod(reference.second+t,24*3600.);
  if (second < 0.) second=second+24*3600.0;
  actual->second=second;

/*   printf("date is: %2d/%2d/%4d \n",actual->day,actual->month,actual->year); */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *sgetnewdate(date_t reference, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  date_t actual;

/* *----------------------------------------------------------------------------
  Get the calendary date corresponding to absolute time t_reference+t */

/*   printf("t is: %f \n",t); */
  getnewdate(reference,t,&actual);

  return(sgetdate(actual));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int clx_initexternal(double time, double dT, state2D_t *set,  mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,k;
  char *datafile=NULL;
  double t,step,cnestime;
  float *buffer[6];
  int mode=LINEAR;
  int nndes=mesh.nvtxs;
  dummy_t memory;

  date_t actual,reference;

  actual= poctime_getdate(time, t_reference, 's', &cnestime);

  exitIfNull(
    datafile=(char *) malloc(strlen(external_path)+strlen(external_root)+64)
    );
  sprintf(datafile,"%s/%s-%4.4d.%2.2d",external_path,external_root,actual.year,actual.month);

  printf("init_external_state: opening %s \n",datafile);

  status=clx_archiveinfo(datafile,&(memory.info));
  if(status!=0) goto error;

  reference=memory.info.reference;
  memory.activated =1;

  status=fe_list(&(memory.info.mesh));
  if(status!=0) goto error;

  status=fe_element_bw(&(memory.info.mesh));
  if(status!=0) goto error;

  status=elements_statev2(memory.info, set, nndes);
  if(status!=0) goto error;

  memory.file=strdup(datafile);

  for(k=0;k<7;k++){
    exitIfNull(
      buffer[k]=(float *) malloc(memory.info.mesh.nvtxs*sizeof(float))
      );
  }
  step=memory.info.sampling;

  switch  (mode) {
    case (LINEAR):
      for (k=0;k<2;k++) {
        exitIfNull(
          set_z[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_u[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_v[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_Hu[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_Hv[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_H[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_sfd[k]=(double *) malloc(nndes*sizeof(double))
          );
      }
/*       t=time; */
      t=time-dT;
      status=clx_readexternal(t,buffer,&set_time[0]);
      status=interpolate_external(set_basic,nndes,set_z[0],buffer[0]);
      status=interpolate_external(set_basic,nndes,set_u[0],buffer[1]);
      status=interpolate_external(set_basic,nndes,set_v[0],buffer[2]);
      status=interpolate_external(set_basic,nndes,set_Hu[0],buffer[3]);
      status=interpolate_external(set_basic,nndes,set_Hv[0],buffer[4]);
      status=interpolate_external(set_basic,nndes,set_H[0],buffer[5]);
      status=interpolate_external(set_basic,nndes,set_sfd[0],buffer[6]);
/*       t=set_time[0]+(double)step; */
      t=time;
      status=clx_readexternal(t,buffer,&set_time[1]);
      status=interpolate_external(set_basic,nndes,set_z[1],buffer[0]);
      status=interpolate_external(set_basic,nndes,set_u[1],buffer[1]);
      status=interpolate_external(set_basic,nndes,set_v[1],buffer[2]);
      status=interpolate_external(set_basic,nndes,set_Hu[1],buffer[3]);
      status=interpolate_external(set_basic,nndes,set_Hv[1],buffer[4]);
      status=interpolate_external(set_basic,nndes,set_H[1],buffer[5]);
      status=interpolate_external(set_basic,nndes,set_sfd[1],buffer[6]);
      break;

    case (QUADRATIC):
      for (k=0;k<3;k++) {
        exitIfNull(
          set_z[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_u[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_v[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_Hu[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_Hv[k]=(double *) malloc(nndes*sizeof(double))
          );
        exitIfNull(
          set_H[k]=(double *) malloc(nndes*sizeof(double))
          );
      }
      t=time-(double)step;
      status=clx_readexternal(t,buffer,&set_time[0]);
      status=interpolate_external(set_basic,nndes,set_z[0],buffer[0]);
      status=interpolate_external(set_basic,nndes,set_u[0],buffer[1]);
      status=interpolate_external(set_basic,nndes,set_v[0],buffer[2]);
      status=interpolate_external(set_basic,nndes,set_Hu[0],buffer[3]);
      status=interpolate_external(set_basic,nndes,set_Hv[0],buffer[4]);
      status=interpolate_external(set_basic,nndes,set_H[0],buffer[5]);
      status=interpolate_external(set_basic,nndes,set_sfd[0],buffer[6]);
      t=set_time[0]+(double)step;
      status=clx_readexternal(t,buffer,&set_time[1]);
      status=interpolate_external(set_basic,nndes,set_z[1],buffer[0]);
      status=interpolate_external(set_basic,nndes,set_u[1],buffer[1]);
      status=interpolate_external(set_basic,nndes,set_v[1],buffer[2]);
      status=interpolate_external(set_basic,nndes,set_Hu[1],buffer[3]);
      status=interpolate_external(set_basic,nndes,set_Hv[1],buffer[4]);
      status=interpolate_external(set_basic,nndes,set_H[1],buffer[5]);
      status=interpolate_external(set_basic,nndes,set_sfd[1],buffer[6]);
      t=set_time[1]+(double)step;
      status=clx_readexternal(t,buffer,&set_time[2]);
      status=interpolate_external(set_basic,nndes,set_z[2],buffer[0]);
      status=interpolate_external(set_basic,nndes,set_u[2],buffer[1]);
      status=interpolate_external(set_basic,nndes,set_v[2],buffer[2]);
      status=interpolate_external(set_basic,nndes,set_Hu[2],buffer[3]);
      status=interpolate_external(set_basic,nndes,set_Hv[2],buffer[4]);
      status=interpolate_external(set_basic,nndes,set_H[2],buffer[5]);
      status=interpolate_external(set_basic,nndes,set_sfd[2],buffer[6]);
      break;
      }

  for(k=0;k<6;k++) free(buffer[k]);

  instantaneous_external_state(time, dT, step, set, nndes);

  if(memory.activated ==1) archive_freeinfo(&(memory.info));
  switch  (mode) {
    case (LINEAR):
      for (k=0;k<2;k++) {
        free(set_z[k]);
        free(set_u[k]);
        free(set_v[k]);
        free(set_Hu[k]);
        free(set_Hv[k]);
        free(set_H[k]);
        }
      break;

    case (QUADRATIC):
      for (k=0;k<3;k++) {
        free(set_z[k]);
        free(set_u[k]);
        free(set_v[k]);
        free(set_Hu[k]);
        free(set_Hv[k]);
        free(set_H[k]);
        }
      break;
      }

  return(0);

 error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int clx_readexternal(double time, float *buffer[5], double *actual_time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
----------------------------------------------------------------------*/
{
  int n,step,status,d1,d2;
  char  *datafile=NULL;
  double start,cnestime;
/*
  meta_archive_t info;
*/
  FILE *in=NULL;
  date_t actual;

  actual= poctime_getdate(time, t_reference, 's', &cnestime);

  exitIfNull(
    datafile=(char *) malloc(strlen(external_path)+strlen(external_root)+64)
    );
  sprintf(datafile,"%s/%s-%4.4d.%2.2d",external_path,external_root,actual.year,actual.month);

  if(obc_memory.file==NULL) {
    obc_memory.file=strdup(datafile);
    status=clx_archiveinfo(datafile,&(obc_memory.info));
    obc_memory.activated =1;
    }

  if(strcmp(datafile,obc_memory.file) !=0) {
    free(obc_memory.file);
    obc_memory.file=strdup(datafile);
    if(obc_memory.activated ==1) archive_freeinfo(&(obc_memory.info));
    status=clx_archiveinfo(datafile,&(obc_memory.info));
/*
    memcpy(&(obc_memory.info),&info,sizeof(meta_archive_t));
*/
    obc_memory.activated =1;
    }
 
  /* model reference */
  d1=julian_day(t_reference);
  /* external archive reference */
  d2=julian_day(obc_memory.info.reference);

  /* shift between model and archive time: archive time=model time -start */
  start=obc_memory.info.start+(d2-d1)*24*3600;

  step=1+ int( (time-start)/obc_memory.info.sampling );
  in=fopen(datafile,"r");
  if(in == NULL) {
    printf("error while opening the datafile %s\n",datafile);
    goto error;
    }

  status=clx_archiveread(in,obc_memory.info,step,buffer,actual_time);
  if(status != 0) {
    printf("error while reading the datafile at step =%d %s\n",datafile,step);
/*     if (step<1) printf("model start before archive content : %s\n", sgetmodeldate(start)); */
    goto error;
    }

  fclose(in);

  /* compute transport */
  for(n=0;n<obc_memory.info.mesh.nvtxs;n++) {
    buffer[3][n]=buffer[1][n]*(obc_memory.info.mesh.vertices[n].h+buffer[0][n]);
    buffer[4][n]=buffer[2][n]*(obc_memory.info.mesh.vertices[n].h+buffer[0][n]);
    buffer[5][n]=obc_memory.info.mesh.vertices[n].h+buffer[0][n];
  /* sfd not here */
    buffer[6][n]=buffer[0][n];
    }

  /* convert archive time in model time */
  *actual_time+=(d2-d1)*24*3600;

  return(status);

error:
  status=-1;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int clx_archiveread(FILE *file,const meta_archive_t & info, int step, float *buffer[3], double *time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  long status;
  long offset;
  int type,nitems,i;
  int nndes,maxn;
  double scale[3],ddum;
  archive1_t *buf1=NULL;
  archive2_t *buf2=NULL;
  archive4_t *buf4=NULL;

  nndes=info.mesh.nvtxs;
  maxn=info.mesh.nnghm;

  /* skip mesh information */

  offset=2*sizeof(int)+nndes*(2*sizeof(float)+maxn*sizeof(int)+sizeof(double));
  status=fseek(file,offset,SEEK_SET);
  if(status !=0) TRAP_ERR_RETURN(errno,1,"fseek(,%d,SEEK_SET) error (%d %s)\n",offset,errno,strerror(errno));

  /* skip time information */

  offset=sizeof(date_t);
  status=fseek(file,offset,SEEK_CUR);
  if(status !=0) TRAP_ERR_RETURN(errno,1,"fseek(,%d,SEEK_CUR) error (%d %s)\n",offset,errno,strerror(errno));

  /* read and skip archive information */

  nitems=fread(&type,sizeof(int),1,file);
  if(nitems!=1) TRAP_ERR_RETURN(errno,1,"fread() returned %d!=1\n",nitems);

#if NEED_SWAP == 1
  type=lswap(type);
#endif

  offset=sizeof(double);
  status=fseek(file,offset,SEEK_CUR);
  if(status !=0) TRAP_ERR_RETURN(errno,1,"fseek(,%d,SEEK_CUR) error (%d %s)\n",offset,errno,strerror(errno));

  /* read scale information */

  nitems=fread(scale,sizeof(double),3,file);

#if NEED_SWAP == 1
  for (i=0;i<3;i++) scale[i]=dswap2(scale[i]);
#endif

  /* skip first time frame */

  if(step > 1) {
    offset=sizeof(double)+nndes*sizeof(archive4_t);
    status=fseek(file,offset,SEEK_CUR);
    if(status !=0) TRAP_ERR_RETURN(errno,1,"fseek(,%d,SEEK_CUR) error (%d %s)\n",offset,errno,strerror(errno));
    }
  else {
    exitIfNull(
      buf4=(archive4_t *) malloc(nndes*sizeof(archive4_t))
      );
    nitems=fread(time,sizeof(double),1,file);
    if(nitems!=1) TRAP_ERR_RETURN(errno,1,"fread() returned %d!=1\n",nitems);
#if NEED_SWAP == 1
    *time=dswap2(*time);
#endif
    nitems=fread(buf4,nndes*sizeof(archive4_t),1,file);
    if(nitems!=1) TRAP_ERR_RETURN(errno,1,"fread() returned %d!=1\n",nitems);
    for (i=0;i<nndes;i++) {
#if NEED_SWAP == 1
/*       buf4[i].z=fswap(buf4[i].z); */
/*       buf4[i].u=fswap(buf4[i].u); */
/*       buf4[i].v=fswap(buf4[i].v); */
      f3swap(&buf4[i].z);
      f3swap(&buf4[i].u);
      f3swap(&buf4[i].v);
#endif
      buffer[0][i]=buf4[i].z;
      buffer[1][i]=buf4[i].u;
      buffer[2][i]=buf4[i].v;
      }
    free(buf4);
    return(0);
    }

  /* skip or read step-2 time frame */

  switch (type) {
    case 1:
      offset=(step-2)*(sizeof(double)+nndes*sizeof(archive1_t));
      status=fseek(file,offset,SEEK_CUR);
      if(status !=0) {
        printf("could not reach variable\n");
        TRAP_ERR_RETURN(errno,1,"fseek(,%d,SEEK_CUR) error (%d %s)\n",offset,errno,strerror(errno));
        }
      break;

    case 2:
      offset=(step-2)*(sizeof(double)+nndes*sizeof(archive2_t));
      status=fseek(file,offset,SEEK_CUR);
      if(status !=0) {
        printf("could not reach variable\n");
        TRAP_ERR_RETURN(errno,1,"fseek(,%d,SEEK_CUR) error (%d %s)\n",offset,errno,strerror(errno));
        }
      break;

    case 4:
      offset=(step-2)*(sizeof(double)+nndes*sizeof(archive4_t));
      status=fseek(file,offset,SEEK_CUR);
      if(status !=0) {
        printf("could not reach variable\n");
        TRAP_ERR_RETURN(errno,1,"fseek(,%d,SEEK_CUR) error (%d %s)\n",offset,errno,strerror(errno));
        }
      break;
    }

  nitems=fread(time,sizeof(double),1,file);
  if(nitems!=1) {
    printf("could not read time\n");
    TRAP_ERR_RETURN(errno,1,"fread() returned %d!=1\n",nitems);
    }
  
#if NEED_SWAP == 1
  *time=dswap2(*time);
#endif

  /* read targeted buffer */

  switch (type) {
    case 1:
      exitIfNull(
        buf1=(archive1_t *) malloc(nndes*sizeof(archive1_t))
        );
      nitems=fread(buf1,nndes*sizeof(archive1_t),1,file);
      if(nitems!=1) {
        printf("could not read variable\n");
        TRAP_ERR_RETURN(errno,1,"fread() returned %d!=1\n",nitems);
        }
      free(buf1);
      break;

    case 2:
      exitIfNull(
        buf2=(archive2_t *) malloc(nndes*sizeof(archive2_t))
        );
      nitems=fread(buf2,nndes*sizeof(archive2_t),1,file);
      if(nitems!=1) {
        printf("could not read variable\n");
        TRAP_ERR_RETURN(errno,1,"fread() returned %d!=1\n",nitems);
        }
      for (i=0;i<nndes;i++) {
#if NEED_SWAP == 1
        buf2[i].z=shswap(buf2[i].z);
        buf2[i].u=shswap(buf2[i].u);
        buf2[i].v=shswap(buf2[i].v);
/*       printf("%d %d \n",i,buf2[i].z); */
#endif
        ddum=buf2[i].z;
        buffer[0][i]=buf2[i].z/scale[0];
        buffer[1][i]=buf2[i].u/scale[1];
        buffer[2][i]=buf2[i].v/scale[2];
/*       printf("%d %d %lf %lf \n",i,buf2[i].u,scale[1],scale[2]); */
        }
      free(buf2);
      break;

    case 4:
      exitIfNull(
        buf4=(archive4_t *) malloc(nndes*sizeof(archive4_t))
        );
      nitems=fread(buf4,nndes*sizeof(archive4_t),1,file);
      if(nitems!=1) {
        printf("could not read variable\n");
        TRAP_ERR_RETURN(errno,1,"fread() returned %d!=1\n",nitems);
        }
      for (i=0;i<nndes;i++) {
#if NEED_SWAP == 1
/*       buf4[i].z=fswap(buf4[i].z); */
/*       buf4[i].u=fswap(buf4[i].u); */
/*       buf4[i].v=fswap(buf4[i].v); */
        f3swap(&buf4[i].z);
        f3swap(&buf4[i].u);
        f3swap(&buf4[i].v);
#endif
        buffer[0][i]=buf4[i].z/scale[0];
        buffer[1][i]=buf4[i].u/scale[1];
        buffer[2][i]=buf4[i].v/scale[2];
        }
      free(buf4);
      break;
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int clx_archivereadheader(FILE *file, meta_archive_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,count,nitems;
  int nndes,maxn;
  int new_format;
  long status;
  char* sdate=NULL;
  date_t date;
  long offset;
  double dum,time;

  float fdum;
  double ddum;
  int idum;
  vertex_t *set=NULL,*seti;
  const range_t<double> degCheck(-1e3,1e3);

  /* read mesh information */
  info->mesh.destroy();

//  info->mesh.degree=1;

  nitems=fread(&idum,sizeof(int),1,file);
  if(nitems!=1) TRAP_ERR_RETURN(-1,1,"fread() returned %d!=1\n",nitems);
#if NEED_SWAP == 1
  idum=lswap (idum);
#endif
  info->mesh.nvtxs=idum;

  nitems=fread(&idum,sizeof(int),1,file);
  if(nitems!=1) TRAP_ERR_RETURN(-1,1,"fread() returned %d!=1\n",nitems);
#if NEED_SWAP == 1
  idum=lswap(idum);
#endif
  info->mesh.nnghm=idum;

  nndes=info->mesh.nvtxs;
  maxn=info->mesh.nnghm;

  exitIfNull(
    info->mesh.vertices=new vertex_t[nndes]
    );
  set= info->mesh.vertices;

  for (i=0; i<nndes; i++) {
    seti=&set[i];
    
    nitems=fread(&fdum,sizeof(float),1,file);
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    seti->lon=fdum*r2d;
    
    nitems=fread(&fdum,sizeof(float),1,file);
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    seti->lat=fdum*r2d;
    
    if( not degCheck.has(seti->lon) or
        not degCheck.has(seti->lat) )
      TRAP_ERR_RETURN(ENOEXEC,1,"node %d has incorrect coordinates (%g,%g)\n",i,seti->lon,seti->lat);
    
    nitems=fread(&ddum,sizeof(double),1,file);
#if NEED_SWAP == 1
    ddum=dswap2(ddum);
#endif
    seti->h=ddum;
    
    exitIfNull(
      seti->ngh=new int[maxn]
      );
    
    for (j=0;j<maxn;j++) {
      nitems=fread(&idum,sizeof(int),1,file);
#if NEED_SWAP == 1
      idum=lswap(idum);
#endif
      seti->ngh[j]=idum;
      }
    }

  new_format=0;
  for (i=0; i<nndes; i++) {
    set[i].nngh=maxn;
    for (j=0;j<maxn;j++) {
      if(set[i].ngh[j]==-1) {
        new_format=1;
        break;
        }
      }
    if(new_format==1) break;
    }

  if(new_format==0) {
    for (i=0; i<nndes; i++) {
      for (j=0;j<maxn;j++) set[i].ngh[j]--;
      }
    }

  for (i=0; i<nndes; i++) {
    set[i].nngh=maxn;
    for (j=0;j<maxn;j++)
      if(set[i].ngh[j]==-1) {
        set[i].nngh=j;
        break;
        }
    if(set[i].nngh==0) {
      printf("trouble with node %d\n",i);
      }
    }

  /* skip time information */

  nitems=fread(&date,sizeof(date_t),1,file);
  if(nitems!=1) TRAP_ERR_RETURN(-1,1,"fread() returned %d!=1\n",nitems);
#if NEED_SWAP == 1
  date=date_swap(date);
#endif
  if(archive_verbose){
    sdate=sgetdate(date);
    printf("archive time reference : %f %s \n",date.second,sdate);
    free(sdate);
    }
/*   t_reference=date; */
  info->reference=date;

  /* read archive information */

  nitems=fread(&idum,sizeof(int),1,file);
  if(nitems!=1) TRAP_ERR_RETURN(-1,1,"fread() returned %d!=1\n",nitems);

#if NEED_SWAP == 1
  idum=lswap(idum);
#endif
/*   archive_type=idum; */
  info->type=idum;

  if((info->type != 1)&&(info->type != 2)&&(info->type != 4))
    TRAP_ERR_RETURN(-1,1,"format error : type is %d, neither 1 nor 2 nor 4\n",info->type);

  nitems=fread(&dum,sizeof(double),1,file);
  if(nitems!=1) TRAP_ERR_RETURN(-1,1,"fread() returned %d!=1\n",nitems);

#if NEED_SWAP == 1
  dum=dswap2(dum);
#endif
  info->sampling=dum;
  /* read scale information */

  nitems=fread(info->scale,sizeof(double),3,file);
  if(nitems!=3) TRAP_ERR_RETURN(-1,1,"fread() returned %d!=3\n",nitems);

#if NEED_SWAP == 1
  for (i=0;i<3;i++) {
    info->scale[i]=dswap2(info->scale[i]);
/*    printf("info->scale[%d]=%lf\n",i,info->scale[i]);*/
    }
#endif

  /* first time frame */

  count=0;
  nitems=fread(&time,sizeof(double),1,file);
  if(nitems!=1)
    TRAP_ERR_RETURN(-1,1,"fread() returned %d!=1 at frame %6d : %12.1lf %s\n",
      nitems,count,time,sgetnewdate(info->reference,time));


#if NEED_SWAP == 1
  time=dswap2(time);
#endif
  info->start=time;

  offset=nndes*sizeof(archive4_t);
  status=fseek(file,offset,SEEK_CUR);
  if(status !=0)
    TRAP_ERR_RETURN(errno,1,"fseek(,%ld,SEEK_CUR) error at frame %6d : %12.1lf %s (%d %s)\n",
      offset,count,time,sgetnewdate(info->reference,time),errno,strerror(errno));
  if(archive_verbose) printf("1st frame %6d : %12.1lf %s \n",count,time,sgetnewdate(info->reference,time));
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int clx_archivereadheader(const char *filename, meta_archive_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  FILE *file;
  
  file=fopen(filename,"r");
  if(file == 0) {
    printf("archive_readheader, failed to open file: %s\n",filename);
    return(-1);
    }
  status=clx_archivereadheader(file, info);
  
  fclose(file);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int clx_archivewriteheader(char *filename, meta_archive_t info, float **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,idum,j,status=0,nitems;
  float fdum;
  FILE *file=NULL;
  date_t date;
  double ddum;

  file=fopen(filename,"w+");

  if(file==NULL) return(-1);

  /* copy the mesh information */

  idum=info.mesh.nvtxs;
#if NEED_SWAP == 1
  idum=lswap(idum);
#endif
  nitems=fwrite(&idum,sizeof(int),1,file);
  if(nitems!=1) return(-1);

  idum=info.mesh.nnghm;
#if NEED_SWAP == 1
  idum=lswap(idum);
#endif
  fwrite(&idum,sizeof(int),1,file);

  for (i=0; i<info.mesh.nvtxs; i++) {
    const vertex_t *vertex=&info.mesh.vertices[i];
    
    fdum=vertex->lon*d2r;
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    
    fdum=vertex->lat*d2r;
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    
    ddum=vertex->h;
#if NEED_SWAP == 1
    ddum=dswap2(ddum);
#endif
    fwrite(&ddum,sizeof(double),1,file);
    
    for (j=0;j<info.mesh.nnghm;j++) {
      if(j<vertex->nngh)
        idum=vertex->ngh[j];
      else
        idum=-1;
#if NEED_SWAP == 1
      idum=lswap(idum);
#endif
      fwrite(&idum,sizeof(int),1,file);
      }
    }

  /* copy time information */

  date=info.reference;
#if NEED_SWAP == 1
  date=date_swap(date);
#endif
  fwrite(&date,sizeof(date_t),1,file);

  /* copy archive information */

  idum=info.type;
#if NEED_SWAP == 1
  idum=lswap(idum);
#endif
  fwrite(&idum,sizeof(int),1,file);

  ddum=info.sampling;
#if NEED_SWAP == 1
  ddum=dswap2(ddum);
#endif
  fwrite(&ddum,sizeof(double),1,file);

  /* copy scales information */

  for (i=0;i<3;i++) {
    ddum=info.scale[i];
#if NEED_SWAP == 1
    ddum=dswap2(ddum);
#endif
    fwrite(&ddum,sizeof(double),1,file);
    }

  /* copy time information */

  ddum=info.start;
#if NEED_SWAP == 1
  ddum=dswap2(ddum);
#endif
  fwrite(&ddum,sizeof(double),1,file);

  for (i=0; i<info.mesh.nvtxs; i++) {
    fdum=buffer[0][i];
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    
    fdum=buffer[1][i];
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    
    fdum=buffer[2][i];
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    }

  fclose(file);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int archive_createheader(char *path,meta_archive_t info,float **buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,idum,j,status=0,nitems;
  float fdum;
  FILE *file=NULL;
  date_t date;
  double ddum;

  file=fopen(path,"w");
  if(file==NULL) return(-1);

  /* copy the mesh information */

  idum=info.mesh.nvtxs;
#if NEED_SWAP == 1
  idum=lswap(idum);
#endif
  nitems=fwrite(&idum,sizeof(int),1,file);
  if(nitems!=1) return(-1);

  idum=info.mesh.nnghm;
#if NEED_SWAP == 1
  idum=lswap(idum);
#endif
  nitems=fwrite(&idum,sizeof(int),1,file);
  if(nitems!=1) return(-1);

  for (i=0; i<info.mesh.nvtxs; i++) {
    fdum=info.mesh.vertices[i].lon*d2r;
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    fdum=info.mesh.vertices[i].lat*d2r;
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    ddum=info.mesh.vertices[i].h;
#if NEED_SWAP == 1
    ddum=dswap2(ddum);
#endif
    fwrite(&ddum,sizeof(double),1,file);
    for (j=0;j<info.mesh.nnghm;j++) {
      idum=info.mesh.vertices[i].ngh[j];
#if NEED_SWAP == 1
      idum=lswap(idum);
#endif
      fwrite(&idum,sizeof(int),1,file);
      }
    }

  /* copy time information */

    date=info.reference;
#if NEED_SWAP == 1
    date=date_swap(date);
#endif
    fwrite(&date,sizeof(date_t),1,file);

  /* copy archive information */

    idum=info.type;
#if NEED_SWAP == 1
    idum=lswap(idum);
#endif
    fwrite(&idum,sizeof(int),1,file);

    ddum=info.sampling;
#if NEED_SWAP == 1
    ddum=dswap2(ddum);
#endif
    fwrite(&ddum,sizeof(double),1,file);

  /* copy scales information */

  for (i=0;i<3;i++) {
    ddum=info.scale[i];
#if NEED_SWAP == 1
    ddum=dswap2(ddum);
#endif
    fwrite(&ddum,sizeof(double),1,file);
    }

  /* copy time information */

  ddum=info.start;
#if NEED_SWAP == 1
  ddum=dswap2(ddum);
#endif
  fwrite(&ddum,sizeof(double),1,file);

  for (i=0; i<info.mesh.nvtxs; i++) {
    fdum=buffer[0][i];
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    fdum=buffer[1][i];
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    fdum=buffer[2][i];
#if NEED_SWAP == 1
    f3swap(&fdum);
#endif
    fwrite(&fdum,sizeof(float),1,file);
    }

  fclose(file);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int archive_skipheader(FILE *file, meta_archive_t info,double t, size_t size)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,count,nitems,found=0,type;
  long status,pos;
  date_t date;
  long offset;
  double dum,time,scale[3];
  int nndes=info.mesh.nvtxs,maxn=info.mesh.nnghm;

  /* skip mesh information */

  offset=2*sizeof(int)+nndes*(2*sizeof(float)+maxn*sizeof(int)+sizeof(double));
  status=fseek(file,offset,SEEK_SET);
  if(status !=0) goto read_error;

  /* skip time information */

  nitems=fread(&date,sizeof(date_t),1,file);
  if(nitems!=1) goto read_error;
#if NEED_SWAP == 1
  date=date_swap(date);
#endif
/*  sdate=sgetdate(date);
  printf("archive time reference : %s \n",sdate);
  free(sdate); */

  nitems=fread(&type,sizeof(int),1,file);
  if(nitems!=1) goto read_error;
#if NEED_SWAP == 1
  type=lswap(type);
#endif
  if(type !=info.type) goto read_error;

  nitems=fread(&dum,sizeof(double),1,file);
  if(nitems!=1) goto read_error;
#if NEED_SWAP == 1
  dum=dswap2(dum);
#endif
  if(dum !=info.sampling) goto read_error;

  /* read scale information */

  nitems=fread(scale,sizeof(double),3,file);
  if(nitems!=3) goto read_error;
#if NEED_SWAP == 1
  for (i=0;i<3;i++) {
    scale[i]=dswap2(scale[i]);
/*    printf("scale[%d]=%lf\n",i,scale[i]);*/
    }
#endif

  /* first time frame */

  count=0;
  nitems=fread(&time,sizeof(double),1,file);
  if(nitems!=1) goto read_error;
#if NEED_SWAP == 1
  time=dswap2(time);
#endif

  if(time==t) {
    printf("need to recreate the archive header: %12.1lf %s \n",time,sgetnewdate(info.reference,time));
    status=1; /*need to recreate the archive header*/
    return(status);
    }

  offset=nndes*sizeof(archive4_t);
  status=fseek(file,offset,SEEK_CUR);
  if(status !=0) goto read_error;
/*   printf("1st frame %6d : %12.1lf %s \n",count,time,sgetnewdate(info.reference,time)); */

  while((!feof(file)) && (found==0)) {
    nitems=fread(&time,sizeof(double),1,file);
    if(nitems==0) break;
#if NEED_SWAP == 1
    time=dswap2(time);
#endif
    count++;
    if((time==t) && (!found)) {
     /*start archiving from that point*/
      pos=ftell(file);
      found=1;
      }
    offset=nndes*size;
    status=fseek(file,offset,SEEK_CUR);
    }

  if(found)  goto update;

/*No previous field for this time frame, check if last found one is AnalysisInterval earlier*/

//   append:
  if(fabs(time+info.sampling-t) > 1.e-3) status=-1;
  else status=0;
/*   printf("append after frame %6d : %12.1lf %s status=%d \n",count,time,sgetnewdate(info.reference,time),status); */

  return(status);

  update:
  status=fseek(file,pos-sizeof(double),SEEK_SET);
  status=2;
  printf("update from frame %6d : %12.1lf %s status=%d \n",count,time,sgetnewdate(info.reference,time),status);
  return(status);

  read_error:
  status=-1;
  printf("read error at frame %6d : %12.1lf %s status=%d \n",count,time,sgetnewdate(info.reference,time),status);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int archive_update2(char *filename, meta_archive_t info, float **buffer, double time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,nitems,status;
  short b;
  FILE *file=NULL;
  archive2_t *archive=NULL;
  date_t date;
  long pos[2];
  
  if(info.type!=2) TRAP_ERR_RETURN(-1,1,"%s coded for type=2, not for type=%d\n",__func__,info.type);
  
  file=fopen(filename,"r+");
  if(file==NULL) TRAP_ERR_RETURN(errno,1,"fopen(\"%s\",\"r+\") error (%d %s)\n",filename,errno,strerror(errno));

  status=archive_skipheader(file,info,time,sizeof(archive2_t));
  if(status!=0) TRAP_ERR_RETURN(status,1,"archive_skipheader((\"%s\"),...) error %d\n",filename,status);

  exitIfNull(
    archive=(archive2_t *) malloc(info.mesh.nvtxs*sizeof(archive2_t))
    );

  for (i=0; i<info.mesh.nvtxs; i++) {
      b=short(NINT(info.scale[0]*buffer[0][i]));
#if NEED_SWAP == 1
      b=shswap(b);
#endif
      archive[i].z=b;
      b=short(NINT(info.scale[1]*buffer[1][i]));
#if NEED_SWAP == 1
      b=shswap(b);
#endif
      archive[i].u=b;
      b=short(NINT(info.scale[2]*buffer[2][i]));
#if NEED_SWAP == 1
      b=shswap(b);
#endif
      archive[i].v=b;
    }
  pos[0]=ftell(file);
/*   printf("position %d, status %d\n",pos[0],status); */
#if NEED_SWAP == 1
  time=dswap2(time);
#endif
  nitems=fwrite(&time,sizeof(double),1,file);
  nitems=fwrite(archive,sizeof(archive2_t)*info.mesh.nvtxs,1,file);
  free(archive);

  fclose(file);
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int clx_archiveinfo(const char *filename, meta_archive_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *file=NULL;  /*  input file */
  int count;
  double time,ddum;
  long status;
  long offset;
  int type,nitems;
  int nndes;
  date_t reference;

  file=fopen(filename,"rb");
  if(file == NULL) TRAP_ERR_RETURN(errno,1,"#error reading file %s (%d %s)\n",filename,errno,strerror(errno));

  status=clx_archivereadheader(file,info);
  if(status!=0) TRAP_ERR_RETURN(status,1,"clx_archivereadheader((\"%s\"),) error %d\n",filename,status);

  count=0;
  nndes=info->mesh.nvtxs;
  type=info->type;

  while(status==0) {
    count++;
    nitems=fread(&ddum,sizeof(double),1,file);
    if(nitems != 1) {
/*       if(archive_verbose) printf("frame %6d : %12.1lf %s\n",count,time,sgetnewdate(info->reference,time),status);  */
      if(archive_verbose) printf("no other time tag after %d (normal file termination)\n",count);
      break;
      }
#if NEED_SWAP == 1
    ddum=dswap2(ddum);
#endif
    time=ddum;
    switch (type) {
      case 1:
        offset=nndes*sizeof(archive1_t);
        break;

      case 2:
        offset=nndes*sizeof(archive2_t);
        break;

      case 4:
        offset=nndes*sizeof(archive4_t);
        break;
      }
    status=fseek(file,offset,SEEK_CUR);
    if(status !=0) {
      printf("truncated time frame at %d (ab-normal file termination)\n",count+1);
      break;
      }
/*      printf("frame %6d : %12.1lf %s\n",count+1,time,sgetnewdate(info->reference,time),status); */
    }

  info->nframe=count;
  fclose(file);
  return(status);
}
