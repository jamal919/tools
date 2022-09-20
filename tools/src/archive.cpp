/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

*******************************************************************************/

#include <config.h>

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"
 
#include "fe.h"
#include "archive.def"
#include "archive.h"
#include "netcdf-proto.h"
#include "poc-time.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int interpolate_external(basic_obsolete_t *basic, int n, double *out, float buffer[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double  z,w;
  int i,j,m;
  
  for(i=0; i<n; i++) {
    z=0;
/*-----------------------------------------------------------------------------
    added to allow tsunami initialisation */
    if(basic[i].element==-1) {
      out[i]=z;
      continue;
      }
    for (j=0;j<3;j++) {
      m=basic[i].nodes[j];
      w=basic[i].beta[j];
      z+=buffer[m]*w;
      }
    out[i]=z;
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int elements_statev1(meta_archive_t info, state2D_t *set, int nndes,mesh_t mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  x,y;
  int k,l,m,n;
  int status;
  int elt;
  int current,nlisted,list[10000];
  date_t reference;
  frame_t frame;
  basic_obsolete_t *set_basic=NULL;

  reference=info.reference;

  exitIfNull(
    set_basic=(basic_obsolete_t *) malloc(nndes*sizeof(basic_obsolete_t))
    );

  status=fe_minmax(info.mesh, frame);

  for(n=0; n<nndes; n++) {
    set_basic[n].element=-1;
    for(k=0;k<3;k++) set_basic[n].nodes[k]=-1;
    for(k=0;k<3;k++) set_basic[n].beta[k]=0.0;
    }

  for(n=0; n<nndes; n++) {
    x=set[n].lon;
    y=set[n].lat;
    if ((x<frame.xmin) || (x>frame.xmax)) continue;
    if ((y<frame.ymin) || (y>frame.ymax)) continue;
/*-----------------------------------------------------------------------------
    initialize element list with element from node neighbours*/
    nlisted=0;
    for(k=0;k<mesh.vertices[n].nngh;k++) {
      m=mesh.vertices[n].ngh[k];
      if(set_basic[m].element!=-1) {
        elt=set_basic[m].element;
        for(l=0;l<nlisted;l++)
          if(list[l]==elt) break;
        list[l]=elt;
        nlisted=l+1;
        }
      }
/*-----------------------------------------------------------------------------
    complete element list with neighbour elements*/
    current=nlisted;
    for(m=0;m<current;m++) {
      for(k=0;k<info.mesh.elt_nngh[list[m]];k++) {
        elt=info.mesh.elt_nghl[list[m]][k];
        for(l=0;l<nlisted;l++)
          if(list[l]==elt) break;
        list[l]=elt;
        updatemax(&nlisted,l+1);
        }
      }
    current=nlisted;
    for(m=0;m<current;m++) {
      for(k=0;k<info.mesh.elt_nngh[list[m]];k++) {
        elt=info.mesh.elt_nghl[list[m]][k];
        for(l=0;l<nlisted;l++)
          if(list[l]==elt) break;
        list[l]=elt;
        updatemax(&nlisted,l+1);
        }
      }

    if(nlisted!=0)
      status=fe_beta_inlist(info.mesh,list, nlisted,x,y,&(set_basic[n].element),set_basic[n].nodes,set_basic[n].beta);
    if(set_basic[n].element==-1)
      status=fe_beta(info.mesh, x, y,&(set_basic[n].element),set_basic[n].nodes,set_basic[n].beta);
/*-----------------------------------------------------------------------------
    commented to allow tsunami initialisation */
/*
    !!! : DESACTIVATED
    if(status!=0) goto error;
*/
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int elements_statev2(meta_archive_t info,state2D_t *set,int nndes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  x,y;
  int k,l,n;
  int status;
  list_t triangles,quadrangles;
  date_t reference;
  grid_t grid;
  basic_obsolete_t *set_basic=NULL;

  reference=info.reference;
  
  exitIfNull(
    set_basic=new basic_obsolete_t [nndes]
    );

  fe_Allocate_and_CreateList(info.mesh,&grid,&triangles);

  for(n=0; n<nndes; n++) {
    set_basic[n].element=-1;
    for(k=0;k<3;k++) set_basic[n].nodes[k]=-1;
    for(k=0;k<3;k++) set_basic[n].beta[k]=0.0;
    }

  for(n=0; n<nndes; n++) {
    x=set[n].lon;
    y=set[n].lat;
    if ((x<grid.xmin) || (x>grid.xmax)) continue;
    if ((y<grid.ymin) || (y>grid.ymax)) continue;
    k=int( floor((x-grid.xmin)/grid.dx) );
    l=int( floor((y-grid.ymin)/grid.dy) );
    status=fe_beta_inlist(info.mesh,triangles.elements[k][l], x,y,&(set_basic[n].element),set_basic[n].nodes,set_basic[n].beta);
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int instantaneous_external_state(double t, double dT, double step, state2D_t *set,int nndes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/*------------------------------------------------------------------------------
  interpolates the external elevation field at open boundary nodes.
  Unit must be cm.
----------------------------------------------------------------------*/

  {
  int i;
  double r;
  double z,u,v,Hu,Hv,H,sfd;

  step=set_time[1]-set_time[0];

  r=(set_time[1]-t)/step;

  for(i=0; i<nndes; i++) {
    z= r*set_z[0][i] +(1.-r)*set_z[1][i];
    u= r*set_u[0][i] +(1.-r)*set_u[1][i];
    v= r*set_v[0][i] +(1.-r)*set_v[1][i];
    Hu=r*set_Hu[0][i]+(1.-r)*set_Hu[1][i];
    Hv=r*set_Hv[0][i]+(1.-r)*set_Hv[1][i];
    H= r*set_H[0][i] +(1.-r)*set_H[1][i];
    sfd= r*set_sfd[0][i] +(1.-r)*set_sfd[1][i];
/*------------------------------------------------------------------------------
    FHL: test to improve external/internal model consistency
    set[i].h=H-z;
    set[i].H=set[i].h+z+sfd;
 */
    set[i].H=set[i].h+sfd;
    set[i].u=u;
    set[i].v=v;
/*     set[i].h=set[i].h-sfd; */
    if(i==130801) {
      printf("%f \n",H);
      }
    }

  r=(set_time[1]-t+dT)/step;

  for(i=0; i<nndes; i++) {
    z= r*set_z[0][i] +(1.-r)*set_z[1][i];
    u= r*set_u[0][i] +(1.-r)*set_u[1][i];
    v= r*set_v[0][i] +(1.-r)*set_v[1][i];
    Hu=r*set_Hu[0][i]+(1.-r)*set_Hu[1][i];
    Hv=r*set_Hv[0][i]+(1.-r)*set_Hv[1][i];
    H= r*set_H[0][i] +(1.-r)*set_H[1][i];
    sfd= r*set_sfd[0][i] +(1.-r)*set_sfd[1][i];
    set[i].Hm=set[i].h+sfd;
    set[i].um=u;
    set[i].vm=v;
    }

  return(0);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void archive_setverbose(int flag)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(flag==0) archive_verbose=0;
  else archive_verbose=1;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archive_read(char *filename, int fmt,const meta_archive_t & info, int step, float *buffer[3], float *mask, double *time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  FILE *file=NULL;
  int h_id,u_id,v_id;
  int verbose=0;
  cdfgbl_t global;
  date_t initial;
  int ncid,time_id;
  size_t start[3],count[3];
  grid_t grid;
  variable_t vinfo;
  
  float mask_[3];
  if(mask==0)
    mask=mask_;
  
  switch (fmt) {
    case TUGO_UG_BINARY_FORMAT:
/* *---------------------------------------------------------------------
      open historical archive files */
      file=fopen(filename,"r");
      if(file==NULL) goto error;
      status=clx_archiveread(file, info, step+1, buffer, time);
      fclose(file);
      break;

    case TUGO_UG_NETCDF_FORMAT:
/* *---------------------------------------------------------------------
      open netcdf UG archive files */
      status=cdf_archiveread(filename, info, step, buffer, time);
      break;

    case TUGO_SG_NETCDF_FORMAT:
/* *---------------------------------------------------------------------
      open netcdf SG archive files */
      status=cdf_globalinfo(filename,&global,verbose);
      status=nc_open(filename,NC_NOWRITE,&ncid);
      status=nc_inq_varid(ncid,"time",&time_id);
      status=nc_inq_varid(ncid,"h",&h_id);
      status=nc_inq_varid(ncid,"u",&u_id);
      status=nc_inq_varid(ncid,"v",&v_id);
      start[0]=step;
      count[0]=1;
      status=nc_get_vara_double(ncid,time_id,start,count,time);
      status=nc_close(ncid);
/*       status=cdf_loadvargrid_2d (filename,h_id,&(info->grid)); */
      status= cdf_loadvar_r1_2d (filename,h_id,0,step,grid,grid.nx,buffer[0],&mask[0],&vinfo);
      status= cdf_loadvar_r1_2d (filename,u_id,0,step,grid,grid.nx,buffer[1],&mask[1],&vinfo);
      status= cdf_loadvar_r1_2d (filename,v_id,0,step,grid,grid.nx,buffer[2],&mask[2],&vinfo);
      break;

    default:
      TRAP_ERR_EXIT(-1,"exiting\n");
      break;
    }
  return(status);

 error:
  status=-1;
  return(status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int archive_freeinfo(meta_archive_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int i;
//   int nndes;
  int status=0;

  if(info==NULL) return(0);

//   nndes=info->mesh.nvtxs;
//
// /*-----------------------------------------------------------------------------
//   exit if structure not allocated */
//   if(info->mesh.vertices==NULL) return(1);
//
//   for (i=0; i<nndes; i++) {
//     /* exit if structure not allocated */
//     if(info->mesh.vertices[i].ngh==NULL) return(-1);
//     free(info->mesh.vertices[i].ngh);
//     }
//
//   free(info->mesh.vertices);
//   info->mesh.vertices=NULL;
//
//   if(info->mesh.triangles!=NULL) free(info->mesh.triangles);
//   info->mesh.triangles=NULL;
  
  info->mesh.destroy();

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int archive_info(const char *filename, int fmt, meta_archive_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  cdfgbl_t global;
  date_t initial;
  int ncid,h_id,time_id,dimid;
  double *time=NULL;
  size_t start[3],count[3];


  switch (fmt) {
    case TUGO_UG_BINARY_FORMAT:
      status=clx_archiveinfo(filename, info);
      break;

    case TUGO_UG_NETCDF_FORMAT:
      status=cdf_archiveinfo(filename, info);
      info->type=2;
      /* NOTE: Oppositely to NetCDF, scale is a DIVIDING parameter!!! */
      aset(info->scale,3,1e3);
      break;

    case TUGO_SG_NETCDF_FORMAT:
      status= cdf_globalinfo(filename,&global,verbose);
      status=nc_open(filename,NC_NOWRITE,&ncid);
      status=nc_inq_varid(ncid,"time",&time_id);
      dimid=global.variable[time_id].dim[0].id;
      exitIfNull(
        time=(double *) malloc(global.dimension[dimid].length*sizeof(double))
        );
      start[0]=0;
      count[0]=global.dimension[dimid].length;
      status=nc_get_vara_double(ncid,time_id,start,count,time);
      initial=poctime_getdatecnes(time[0],'s');
      info->reference=poctime_getdatecnes(0.0,'s');
      info->nframe=global.dimension[dimid].length;
      info->sampling=time[1]-time[0];
      info->start=time[0];
      status=nc_inq_varid(ncid,"h",&h_id);
      status=nc_close(ncid);
/*       status=cdf_loadvargrid_2d (filename,h_id,&(info->grid)); */
      status=poc_getgrid2d (filename,global, global.variable[h_id],&(info->grid));
      break;

    default:
      TRAP_ERR_EXIT(-1,"exiting\n");
      break;
    }
  return(status);
}

/*----------------------------------------------------------------------------*/

  int fe_framemapr1(char *input, meta_archive_t info, grid_t grid, float mask, int *elts, int step, float *buf, double *tag)

/*----------------------------------------------------------------------------*/
{
  FILE *in=NULL;
  float *buffer[3]={NULL,NULL,NULL};
  int status;
  double time;
  int k,nndes=info.mesh.nvtxs;

  in=fopen(input, "rb");
  if(in==0) TRAP_ERR_EXIT(-1,"file opening issue : %s \n",input);


  for(k=0;k<3;k++) {
    buffer[k]=(float *)malloc(nndes*sizeof(float));
    if(buffer[k] == NULL) {
      printf("#fe_framemapr1 memory allocation error for buffer N= %d \n",nndes);
      goto error;
      }
    }

  status=clx_archiveread(in, info, step, buffer, &time);
  if(status != 0) {
    printf("#fe_framemapr1 read error at frame= %d \n",step);
    goto error;
    }
  printf("#treating frame for %s\n",sgetnewdate(info.reference,time));

  status=fe_map(info.mesh,buffer[0],grid,elts,buf,mask);
  *tag=time;

  for(k=0;k<3;k++) if(buffer[k] != NULL) free(buffer[k]);
  fclose(in);
  return(0);

error:
  for(k=0;k<3;k++) if(buffer[k] != NULL) free(buffer[k]);
  if(in != NULL)  fclose(in);
  return(-1);
}

/*----------------------------------------------------------------------------*/

  int fe_framemapr2(char *input, meta_archive_t info, grid_t grid, float mask, int *elts, int step, float *bufx, float *bufy, double *tag)

/*----------------------------------------------------------------------------*/
{
  FILE *in=NULL;  			/*  input file */
  float *buffer[3]={NULL,NULL,NULL};
  int status;
  double time;
  int k,nndes=info.mesh.nvtxs;
  
  in=fopen(input, "rb");
  if(in==0) TRAP_ERR_EXIT(-1,"file opening issue : %s \n",input);
  
  for(k=0;k<3;k++) {
    buffer[k]=(float *)malloc(nndes*sizeof(float));
    if(buffer[k] == NULL) {
      printf("#model2d_map1d memory allocation error for buffer N= %d \n",nndes);
      goto error;
      }
    }

  status=clx_archiveread(in, info, step, buffer, &time);
  if(status != 0) {
    printf("#model2d_map1d read error at frame= %d \n",step);
    goto error;
    }
  printf("#treating frame for %s\n",sgetnewdate(info.reference,time));

  status=fe_map(info.mesh,buffer[1],grid,elts,bufx,mask);
  status=fe_map(info.mesh,buffer[2],grid,elts,bufy,mask);
  *tag=time;

  for(k=0;k<3;k++) if(buffer[k] != NULL) free(buffer[k]);
  fclose(in);
  return(0);

error:
  for(k=0;k<3;k++) if(buffer[k] != NULL) free(buffer[k]);
  if(in != NULL)  fclose(in);
  return(-1);
}

