#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
 
#include "fe.h"
#include "archive.def"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_archiveinfo(const char *filename, meta_archive_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,dim,time_id;
  int verbose=0;
  double time[2];
  size_t start[10],count[10];
  int ncid;

  cdfgbl_t g_info;
//   cdfvar_t v_info;

  status= cdf_globalinfo(filename,&g_info,verbose);
  
  if(status!=0) TRAP_ERR_EXIT(-1,"getting netcdf information on file=%s failed\n", filename);
  
  time_id=cdf_identify(g_info, (char *)("time") );
  
  status=nc_open(filename,NC_NOWRITE,&ncid);
  
  status= poc_gettime (ncid, &info->reference, 0, &info->nframe);
  NC_TRAP_ERROR(wexit,status,verbose>=0,"error getting time origin in file=%s", filename);
  
  start[0]=0;
  count[0]=2;
  status=nc_get_vara_double(ncid,time_id,start,count,time);

  status=nc_close(ncid);

  info->sampling=time[1]-time[0];
  info->start=time[0];

/*   status= cdf_varinfo(filename,vdata, &v_info);  */
  
  status=fe_readmesh3d(filename,&(info->mesh),verbose);
  NC_TRAP_ERROR(return,status,verbose>0,"fe_readmesh3d(\"%s\",,) error", filename);

  for(dim=0;dim<g_info.ndimsp;dim++) {
    if(isT(g_info.dimension[dim].name)) {
      info->nframe=g_info.dimension[dim].length;
      continue;
      }
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_archiveread(const char *filename,const meta_archive_t & info, int frame, float **buffer, double *time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  date_t actual;
  int ncid;
  int time_id;
  int h_id,u_id,v_id;
  const char *h_name=0,*u_name=0,*v_name=0;
  size_t start[10],count[10];

  status=nc_open(filename,NC_NOWRITE,&ncid);
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_open error");
  
/*------------------------------------------------------------------------------
  collect time */
  
  status=nc_inq_varid(ncid,"time",&time_id);
  if(status!=0) time_id=-1;

  start[0]=frame;
  count[0]=1;
  status=nc_get_vara_double(ncid,time_id,start,count,time);
  
/*------------------------------------------------------------------------------
  set variables' names */
  
  switch(info.flag){
    case FLAG_STATE:
      h_name="elevation";
      u_name="ubar";
      v_name="vbar";
      break;
  
    case FLAG_FORCING:
      h_name="Pa";
      status=nc_inq_varid(ncid,h_name,&h_id);
      if(status!=0) {
        h_name="ibd";
        }
      u_name="wsx";
      v_name="wsy";
      break;
  
    default:
      TRAP_ERR_EXIT(ENOEXEC,"not coded for flag=%d\n",info.flag);
    }
  
/*------------------------------------------------------------------------------
  collect variables' IDs */
  
  status=nc_inq_varid(ncid,h_name,&h_id);
  if(status!=0) h_id=-1;
  
  status=nc_inq_varid(ncid,u_name,&u_id);
  if(status!=0) u_id=-1;
  
  status=nc_inq_varid(ncid,v_name,&v_id);
  if(status!=0) v_id=-1;
  
  status=nc_close(ncid);
  
/*------------------------------------------------------------------------------
  read data */
  
//  printf("load %s, frame=%d\n",filename,frame);
  
  if(h_id>=0){
    status=poc_get_UG3D(filename,frame,h_name,buffer[0],0);
    NC_TRAP_ERROR(return,status,1,"poc_get_UG3D error");
    }
  else{
    aset(buffer[0],info.mesh.nvtxs,0.f);
    }
  
  if(u_id>=0){
    status=poc_get_UG3D(filename,frame,u_name,buffer[1],0);
    NC_TRAP_ERROR(return,status,1,"poc_get_UG3D error");
    }
  else{
    aset(buffer[1],info.mesh.nvtxs,0.f);
    }
  
  if(v_id>=0){
    status=poc_get_UG3D(filename,frame,v_name,buffer[2],0);
    NC_TRAP_ERROR(return,status,1,"poc_get_UG3D error");
    }
  else{
    aset(buffer[2],info.mesh.nvtxs,0.f);
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_readexternal(double time, date_t reference, float *buffer[5], double *actual_time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,step,status,d1,d2;
  char  *datafile=NULL;
  double first,cnestime;
  date_t actual;
  int ncid;
  int time_id;
  int h_id;
  int u_id;
  int v_id;
  int p_id;
  int sfd_id;
  size_t start[10],count[10];
  float mask;

  actual= poctime_getdate(time, reference, 's', &cnestime);

  exitIfNull(
    datafile=(char *) malloc(strlen(external_path)+strlen(external_root)+64)
    );
  sprintf(datafile,"%s/%s-%4.4d.%2.2d.nc",external_path,external_root,actual.year,actual.month);

  if(obc_memory.file==NULL) {
    obc_memory.file=strdup(datafile);
    status=cdf_archiveinfo(datafile,&(obc_memory.info));
    obc_memory.activated =1;
    }

  if(strcmp(datafile,obc_memory.file) !=0) {
    free(obc_memory.file);
    obc_memory.file=strdup(datafile);
    if(obc_memory.activated ==1) archive_freeinfo(&(obc_memory.info));
    status=cdf_archiveinfo(datafile,&(obc_memory.info));
    obc_memory.activated =1;
    }
 
/*------------------------------------------------------------------------------
  extraction reference */
  d1=julian_day(reference);
/*------------------------------------------------------------------------------
  external archive reference */
  d2=julian_day(obc_memory.info.reference);
/*------------------------------------------------------------------------------
  shift between extraction and archive time: archive time=model time -start */
  first=obc_memory.info.start+(d2-d1)*24*3600;

  step=int( floor((time-first)/obc_memory.info.sampling+0.5) );

  status=nc_open(datafile,NC_NOWRITE,&ncid);
  if(status!=0) goto error;

  status=nc_inq_varid(ncid,"time",&time_id);
  if(status!=0) time_id=-1;

  start[0]=step;
  count[0]=1;
  status=nc_get_vara_double(ncid,time_id,start,count,actual_time);
/*------------------------------------------------------------------------------
  convert archive time in extraction time */
  *actual_time+=(d2-d1)*24*3600.;

  status=nc_inq_varid(ncid,"elevation",&h_id);
  if(status!=0) h_id=-1;

  status=nc_inq_varid(ncid,"ubar",&u_id);
  if(status!=0) u_id=-1;

  status=nc_inq_varid(ncid,"vbar",&v_id);
  if(status!=0) v_id=-1;

  status=nc_inq_varid(ncid,"pmsl",&p_id);
  if(status!=0) p_id=-1;

  status=nc_inq_varid(ncid,"sfd",&sfd_id);
  if(status!=0) sfd_id=-1;

  status=nc_close(ncid);

  printf("load %s/%s-%4.4d.%2.2d.nc, frame=%d\n",external_path,external_root,actual.year,actual.month,step);

/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=poc_getshort_nt(datafile,obc_memory.info.mesh,step,h_id,buffer[0],&mask);
  if(status!=0) goto error;
/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=poc_getshort_nt(datafile,obc_memory.info.mesh,step,u_id,buffer[1],&mask);
  if(status!=0) goto error;
/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=poc_getshort_nt(datafile,obc_memory.info.mesh,step,v_id,buffer[2],&mask);
  if(status!=0) goto error;

  if(sfd_id!=-1) {
/**----------------------------------------------------------------------------
    obsolete call, will be suppressed */
    status=poc_getshort_nt(datafile,obc_memory.info.mesh,step,sfd_id,buffer[6],&mask);
    if(status!=0) goto error;
    }
  else
    for(n=0;n<obc_memory.info.mesh.nvtxs;n++) buffer[6][n]=buffer[0][n];

  for(n=0;n<obc_memory.info.mesh.nvtxs;n++) {
    buffer[3][n]=buffer[1][n]*(obc_memory.info.mesh.vertices[n].h+buffer[0][n]);
    buffer[4][n]=buffer[2][n]*(obc_memory.info.mesh.vertices[n].h+buffer[0][n]);
    buffer[5][n]=obc_memory.info.mesh.vertices[n].h+buffer[0][n];
    }
  
  return(status);

 error:
  status=-1;
  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_initexternal(double time, double dT, state2D_t *set,  mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,k;
  char  *datafile=NULL;
  double t,step,cnestime,first;
  float *buffer[10];
  int mode=LINEAR;
  int nndes=mesh.nvtxs;
  int ncid;
  int time_id;
  int h_id;
  int u_id;
  int v_id;
  int p_id;
  int sfd_id;
  size_t start[10],count[10];

  date_t actual,reference;

  actual= poctime_getdate(time, t_reference, 's', &cnestime);

  exitIfNull(
    datafile=(char *) malloc(strlen(external_path)+strlen(external_root)+64)
    );
  sprintf(datafile,"%s/%s-%4.4d.%2.2d.nc",external_path,external_root,actual.year,actual.month);

  printf("init_external_state: opening %s \n",datafile);

  status=cdf_archiveinfo(datafile,&(obc_memory.info));
  if(status!=0) goto error;

  reference=obc_memory.info.reference;
  obc_memory.activated =1;

  status=fe_list(&(obc_memory.info.mesh));
  if(status!=0) goto error;

  status=fe_element_bw(&(obc_memory.info.mesh));
  if(status!=0) goto error;

  status=elements_statev2(obc_memory.info, set, nndes);
  if(status!=0) goto error;

  obc_memory.file=strdup(datafile);

  for(k=0;k<7;k++) {
    exitIfNull(
      buffer[k]=(float *) malloc(obc_memory.info.mesh.nvtxs*sizeof(float))
      );
  }
  step=obc_memory.info.sampling;

  status=nc_open(datafile,NC_NOWRITE,&ncid);

  status=nc_inq_varid(ncid,"time",&time_id);
  if(status!=0) time_id=-1;

  start[0]=0;
  count[0]=1;
  status=nc_get_vara_double(ncid,time_id,start,count,&first);
/*   if(status==0) initial=time_getcnesdate(first,'s'); */

  status=nc_inq_varid(ncid,"elevation",&h_id);
  if(status!=0) h_id=-1;

  status=nc_inq_varid(ncid,"ubar",&u_id);
  if(status!=0) u_id=-1;

  status=nc_inq_varid(ncid,"vbar",&v_id);
  if(status!=0) v_id=-1;

  status=nc_inq_varid(ncid,"pmsl",&p_id);
  if(status!=0) p_id=-1;

  status=nc_inq_varid(ncid,"sfd",&sfd_id);
  if(status!=0) sfd_id=-1;

  status=nc_close(ncid);

  switch  (mode)
    {
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
      status=cdf_readexternal(t,t_reference,buffer,&set_time[0]);
      status=interpolate_external(set_basic,nndes,set_z[0],buffer[0]);
      status=interpolate_external(set_basic,nndes,set_u[0],buffer[1]);
      status=interpolate_external(set_basic,nndes,set_v[0],buffer[2]);
      status=interpolate_external(set_basic,nndes,set_Hu[0],buffer[3]);
      status=interpolate_external(set_basic,nndes,set_Hv[0],buffer[4]);
      status=interpolate_external(set_basic,nndes,set_H[0],buffer[5]);
      status=interpolate_external(set_basic,nndes,set_sfd[0],buffer[6]);
/*       t=set_time[0]+(double)step; */
      t=time-dT+obc_memory.info.sampling;
      status=cdf_readexternal(t,t_reference,buffer,&set_time[1]);
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
      status=cdf_readexternal(t,t_reference,buffer,&set_time[0]);
      status=interpolate_external(set_basic,nndes,set_z[0],buffer[0]);
      status=interpolate_external(set_basic,nndes,set_u[0],buffer[1]);
      status=interpolate_external(set_basic,nndes,set_v[0],buffer[2]);
      status=interpolate_external(set_basic,nndes,set_Hu[0],buffer[3]);
      status=interpolate_external(set_basic,nndes,set_Hv[0],buffer[4]);
      status=interpolate_external(set_basic,nndes,set_H[0],buffer[5]);
      status=interpolate_external(set_basic,nndes,set_sfd[0],buffer[6]);
      t=set_time[0]+(double)step;
      status=cdf_readexternal(t,t_reference,buffer,&set_time[1]);
      status=interpolate_external(set_basic,nndes,set_z[1],buffer[0]);
      status=interpolate_external(set_basic,nndes,set_u[1],buffer[1]);
      status=interpolate_external(set_basic,nndes,set_v[1],buffer[2]);
      status=interpolate_external(set_basic,nndes,set_Hu[1],buffer[3]);
      status=interpolate_external(set_basic,nndes,set_Hv[1],buffer[4]);
      status=interpolate_external(set_basic,nndes,set_H[1],buffer[5]);
      status=interpolate_external(set_basic,nndes,set_sfd[1],buffer[6]);
      t=set_time[1]+(double)step;
      status=cdf_readexternal(t,t_reference,buffer,&set_time[2]);
      status=interpolate_external(set_basic,nndes,set_z[2],buffer[0]);
      status=interpolate_external(set_basic,nndes,set_u[2],buffer[1]);
      status=interpolate_external(set_basic,nndes,set_v[2],buffer[2]);
      status=interpolate_external(set_basic,nndes,set_Hu[2],buffer[3]);
      status=interpolate_external(set_basic,nndes,set_Hv[2],buffer[4]);
      status=interpolate_external(set_basic,nndes,set_H[2],buffer[5]);
      status=interpolate_external(set_basic,nndes,set_sfd[2],buffer[6]);
      break;
      }

  for(k=0;k<7;k++) free(buffer[k]);

  instantaneous_external_state(time, dT, step, set, nndes);

  if(obc_memory.activated ==1) archive_freeinfo(&(obc_memory.info));
  switch  (mode)
    {
    case (LINEAR):
      for (k=0;k<2;k++) {
        free(set_z[k]);
        free(set_u[k]);
        free(set_v[k]);
        free(set_Hu[k]);
        free(set_Hv[k]);
        free(set_H[k]);
        free(set_sfd[k]);
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
        free(set_sfd[k]);
        }
      break;
      }


  return(0);

 error:
  return(-1);

}
