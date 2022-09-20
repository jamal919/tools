

/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <sys/timeb.h>

#if NEED_IEEEFP == 1
#include <ieeefp.h>
#endif

#include "tools-structures.h"

#include "netcdf-proto.h"
#include "poc-netcdf-data.hpp"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_inquire_format (const char* input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int verbose;
  cdfgbl_t grid_info;
  int  fmt;

  verbose=0;
  status=cdf_globalinfo(input,&grid_info,verbose);
  if(status!=0) return(-1);

  fmt=0;  /** old format */

  for(n=0;n<grid_info.ngattsp;n++) {
    if(strcmp(grid_info.attribute[n].name,"Conventions")==0) {
      fmt=1;  /** new format */
      break;
      }
    }
  
  grid_info.destroy();

  return(fmt);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_loadgrid_01 (const char* filename,grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, ncid;
  size_t index[1];
/*   cdf_fheader header; */
  int ndimsp,nvarsp,ngattsp,unlimdimidp,var,nattsp,dimids[50];
  nc_type xtypep;
  char **name=NULL;

  status=nc_open(filename,0,&ncid);
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_open(\"%s\",0,) error",filename);

  status=nc_inq(ncid,&ndimsp,&nvarsp,&ngattsp,&unlimdimidp);
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_inq((\"%s\"),...) error",filename);

  exitIfNull(
    name=new char*[nvarsp]
    );
  for (var=0;var<nvarsp;var++) {
    exitIfNull(
      name[var]=new char[1024]
      );
    status=nc_inq_var(ncid,var,name[var],&xtypep,&ndimsp,dimids,&nattsp);
    if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_inq_var((\"%s\"),%d,...) error",filename,var);
    }

  if (strcmp(name[0],"x_range")!=0) TRAP_ERR_RETURN(-1,"first variable is %s, not x_range.\n",name[0]);
  
  deletep2D(&name,nvarsp);

  var=0;
  index[0]=0;
  status=nc_get_var1_double(ncid,var,index,&(grid->xmin));
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_get_var1_double((\"%s\"),%d,...) error",filename,var);
  index[0]=1;
  status=nc_get_var1_double(ncid,var,index,&(grid->xmax));
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_get_var1_double((\"%s\"),%d,...) error",filename,var);

  var=1;
  index[0]=0;
  status=nc_get_var1_double(ncid,var,index,&(grid->ymin));
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_get_var1_double((\"%s\"),%d,...) error",filename,var);
  index[0]=1;
  status=nc_get_var1_double(ncid,var,index,&(grid->ymax));
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_get_var1_double((\"%s\"),%d,...) error",filename,var);

  var=3;
  index[0]=0;
  status=nc_get_var1_double(ncid,var,index,&(grid->dx));
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_get_var1_double((\"%s\"),%d,...) error",filename,var);
  index[0]=1;
  status=nc_get_var1_double(ncid,var,index,&(grid->dy));
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_get_var1_double((\"%s\"),%d,...) error",filename,var);

  var=4;
  index[0]=0;
  status=nc_get_var1_int(ncid,var,index,&(grid->nx));
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_get_var1_int((\"%s\"),%d,...) error",filename,var);
  index[0]=1;
  status=nc_get_var1_int(ncid,var,index,&(grid->ny));
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_get_var1_int((\"%s\"),%d,...) error",filename,var);

  grid->modeH=0;
  status=nc_close(ncid);
  if(status!=0) NC_TRAP_ERROR(return,status,1,"nc_close((\"%s\")) error",filename);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_loadgrid_02 (const char* input,grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
int    /*v,*/status;
int    verbose;
cdfgbl_t data_info,grid_info;
cdfvar_t h_info;

  verbose=0;
  status=cdf_globalinfo(input,&data_info,verbose);
  if(status != 0) return(-1);
//   for (v=0;v<data_info.nvarsp;v++) {
//     printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
//     }
  status=cdf_varinfo(input,"z",&h_info,1);
  if(status != 0) return(-1);

  status=cdf_globalinfo(input,&grid_info,verbose);
  if(status != 0) return(-1);
  
  status= poc_getgrid2d (input, grid_info, h_info, grid);
  if(status != 0) return(-1);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_loadgrid (const char* input,grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;
  int  fmt;

  fmt=grd_inquire_format (input);

/*------------------------------------------------------------------------------
  distinguish old and new GMT format */
  switch (fmt) {
    case 0:
      status=grd_loadgrid_01(input,grid);
      break;

    case 1:
//      status=grd_loadgrid_02(input,grid);
      status= poc_get_grid(input,"z",grid,0);
      break;

    default:
      return(-1);
      break;
    }
 return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_savegrid(const char* filename, grid_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, ncid;
  size_t lengthhp[2],index[10];
  int ndimsp,nvarsp,ngattsp,unlimdimidp,dim,var,att,length,dimids[50];
  char *s1=(char *)("side"),*s2=(char *)("xysize");
  char *varname[6],*dimname[2], *attname[2], *attvalue[2];
  nc_type vartype[6];
  int varid[6],vardim[6],varndimsp[6];
  int i;
  double ddum;

//   status=nc_create(filename,NC_CLOBBER,&ncid);
//   status=nc_create(filename,NC_CDF5,&ncid); 
  status=nc_create(filename,NC_NETCDF4,&ncid); 
  if(status != NC_NOERR) goto error;
  
  ndimsp=2;

  dimname[0]=strdup(s1);
  lengthhp[0]=2;

  dimname[1]=strdup(s2);
  lengthhp[1]=(size_t) grid.nx*(size_t) grid.ny;

  nvarsp=5; /* z not defined here */
  ngattsp=2;
  unlimdimidp=-1;
 
  for (dim=0;dim<ndimsp;dim++) {
    status=nc_def_dim(ncid,dimname[dim],lengthhp[dim],&dimids[dim]);
    if(status != NC_NOERR) goto error;
    }

  varname[0]=(char *)("x_range");
  vartype[0]=NC_DOUBLE;
  varndimsp[0]= 1;
  vardim[0]= 0;
  varname[1]=(char *)("y_range");
  vartype[1]=NC_DOUBLE;
  varndimsp[1]= 1;
  vardim[1]= 0;
  varname[2]=(char *)("z_range");
  vartype[2]=NC_DOUBLE;
  varndimsp[2]= 1;
  vardim[2]= 0;
  varname[3]=(char *)("spacing");
  vartype[3]=NC_DOUBLE;
  varndimsp[3]= 1;
  vardim[3]= 0;
  varname[4]=(char *)("dimension");
  vartype[4]=NC_INT;
  varndimsp[4]= 1;
  vardim[4]= 0;
/*
  varname[5]="z";
  vartype[5]=NC_FLOAT;
  varndimsp[5]= 1;
  vardim[5]= 1;
*/
  for (var=0;var<nvarsp;var++) {
    status=nc_def_var(ncid,varname[var],vartype[var],varndimsp[var],&vardim[var],&varid[var]);
    if(status != NC_NOERR) goto error;
    }

  ngattsp=2;

  att=0;
  attname[0]=strdup("title");
  attvalue[0]=strdup("\0");
  length=strlen(attvalue[att])+1;
  length=1;
  status=nc_put_att_text(ncid,NC_GLOBAL,attname[att],length,attvalue[att]);
  if(status != NC_NOERR) goto error;

  att=1;
  attname[1]=strdup("source");
  exitIfNull(
    attvalue[1]=(char *) malloc(481)
    );
  for (i=0;i<481;i++) *(attvalue[1]+i)=0;
  strcpy(attvalue[1],"created by grd tools library, F. Lyard 2008");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,NC_GLOBAL,attname[att],length,attvalue[att]);

  time_t creation_time;
  status=time(&creation_time);

  attname[1]=strdup("creation");
  for (i=0;i<481;i++) *(attvalue[1]+i)=0;
  strcpy(attvalue[1],ctime(&creation_time));
  length=strlen(attvalue[att]);
  *(attvalue[1]+length-1)=0;
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,NC_GLOBAL,attname[att],length,attvalue[att]);


  att=0;
  attname[att]=strdup("units");
  attvalue[att]=strdup("user_x_unit");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,0,attname[att],length,attvalue[att]);
  if(status != NC_NOERR) goto error;

  free(attvalue[att]);
  attvalue[att]=strdup("user_y_unit");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,1,attname[att],length,attvalue[att]);
  if(status != NC_NOERR) goto error;

  free(attvalue[att]);
  attvalue[att]=strdup("user_z_unit");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,2,attname[att],length,attvalue[att]);
  if(status != NC_NOERR) goto error;

  status=nc_enddef(ncid);
  if(status != NC_NOERR) goto error;

  var=0;
  index[0]=0;
  ddum=grid.xmin;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  if(status != NC_NOERR) goto error;
  index[0]=1;
  ddum=grid.xmax;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  if(status != NC_NOERR) goto error;

  var=1;
  index[0]=0;
  ddum=grid.ymin;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  if(status != NC_NOERR) goto error;
  index[0]=1;
  ddum=grid.ymax;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  if(status != NC_NOERR) goto error;

  var=2;
  index[0]=0;
  ddum=0;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  if(status != NC_NOERR) goto error;
  index[0]=1;
  ddum=1;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  if(status != NC_NOERR) goto error;

  var=3;
  index[0]=0;
  ddum=grid.dx;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  if(status != NC_NOERR) goto error;
  index[0]=1;
  ddum=grid.dy;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  if(status != NC_NOERR) goto error;

  var=4;
  index[0]=0;
  status=nc_put_var1_int(ncid,var,index,&(grid.nx));
  if(status != NC_NOERR) goto error;
  index[0]=1;
  status=nc_put_var1_int(ncid,var,index,&(grid.ny));
  if(status != NC_NOERR) goto error;

  status=nc_close(ncid);
  if(status != NC_NOERR) goto error;

  return(status);

  error:
  if(status != NC_NOERR) {
    nc_advise("ncattget", status, "ncid %d", ncid);
    return(status);
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void grd_getdimension(const char *filename,int* nv,int* ni,int* nj,int* nk,int* nd,int* nt,int* status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  grid_t grid;

  *status=grd_loadgrid (filename, &grid);
  if(*status != 0) return;
  *nv=1;
  *nk=1;
  *nt=1;
  *nd=1;

  *ni=grid.nx;
  *nj=grid.ny;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_checkfile (const char *filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  grid_t grid;
  int status;

  status=grd_loadgrid (filename, &grid);

  if(status == 0) {
    printf("grd_checkfile--------------------------------------  \n");
    status=cdf_info(filename);
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_DuplicateTag(const char* in, const char* out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, cerr, ncid;
  int varid;
  poc_var_t var;
  short *buf=0;
  size_t size=1;
  size_t count;

  status=nc_open(in,0,&ncid);
  if (status !=NC_NOERR) {
    status=-1;
    return(status);
    }
  
  status=poc_inq_var(ncid,"tag",&var,0);
  if(status!=0) goto terminate;
    
  for(int k=0;k<var.dimensions.size();k++) size*=var.dimensions[k].len;
  
  buf=new short[size];
  status=poc_get_var(ncid,var.id,buf);
  if(status!=0) goto terminate;
   
  count=occurence<short>((short) -1, buf, size);
  status=nc_close(ncid);

  status=nc_open(out,NC_WRITE,&ncid);
  if (status !=NC_NOERR) {
    delete[] buf;
    status=-1;
    return(status);
    }
  status=poc_def_var(ncid,var,&varid);
  if(status!=0) goto terminate;

  status=poc_put_var(ncid, varid, buf);
  if(status!=0) goto terminate;
  
  delete[] buf;
  
  printf("tag detected and copied\n");

terminate:

  cerr=nc_close(ncid);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int grd_loadr1_01_template(const char* filename, grid_t grid, T *buf, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  handle new grd format, no posterior mirroring required

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status, ncid;
  int varid;

  status=nc_open(filename,0,&ncid);
  if (status !=NC_NOERR) {
    status=-1;
    return(status);
    }

  varid=5;
  status=poc_get_var(ncid,varid,buf);

  status = poc_get_att(ncid, varid, "_FillValue", mask);
  if(status != NC_NOERR)
    *mask = 1.0e+10;

  status=nc_close(ncid);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_loadr1_01 (const char* filename, grid_t grid, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=grd_loadr1_01_template(filename, grid, buf, mask);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_loadr1_01 (const char* filename, grid_t grid, short *buf, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=grd_loadr1_01_template(filename, grid, buf, mask);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int grd_loadr1_02_template(const char* input, grid_t grid, T *buf, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  handle new grd format, no posterior mirroring required

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int    verbose,frame=0;
  cdfgbl_t data_info,grid_info;
  cdfvar_t h_info;

  verbose=0;

  status=cdf_globalinfo(input,&data_info,verbose);
  status=cdf_varinfo(input,"z",&h_info,1);

  status= poc_getvar2d(input, h_info.id, frame, buf, mask, h_info);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int grd_loadr1_02 (const char* input, grid_t grid, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=grd_loadr1_02_template(input, grid, buf, mask);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int grd_loadr1_02 (const char* input, grid_t grid, short *buf, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=grd_loadr1_02_template(input, grid, buf, mask);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int grd_loadr1_template(const char* input, grid_t grid, T *buf, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  fmt,status;

  fmt=grd_inquire_format (input);

  switch (fmt) {
    case 0:
      status=grd_loadr1_01(input, grid, buf, mask);
      status= grd_mirror_r(grid, grid.nx, buf, *mask);
      break;

    case 1:
      status=grd_loadr1_02(input, grid, buf, mask);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_loadr1 (const char* input, grid_t grid, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  status;

  status=grd_loadr1_template(input, grid, buf, mask);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_loadr1 (const char* input, grid_t grid, short *buf, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  status;

  status=grd_loadr1_template(input, grid, buf, mask);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_loads1 (const char* filename, grid_t grid, int n, short *buf, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, ncid, varid;

  status=nc_open(filename,0,&ncid);
  if (status !=NC_NOERR) {
    status=-1;
    return(status);
    }

  varid=5;
  status=nc_get_var_short(ncid,varid,buf);

  status = nc_get_att(ncid, varid, "_FillValue", mask);
  if(status != NC_NOERR)
    *mask = 256*127+255;


  status=nc_close(ncid);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int grd_mirror_r_template (grid_t grid, int n, T *buffer, T mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T tmp;
  size_t i,j,m1,m2;

  for (j=0;j<grid.ny/2;j++) {
    for (i=0;i<grid.nx;i++) {
      m1=grid.nx*j+i;
      m2=grid.nx*(grid.ny-j-1)+i;
      tmp=buffer[m1];
      buffer[m1]=buffer[m2];
      buffer[m2]=tmp;
/*
      if(isnanf(buffer[m1])) buffer[m1]=mask;
      if(isnanf(buffer[m2])) buffer[m2]=mask;
*/
      }
    }
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_mirror_r (grid_t grid, int n, float *buffer, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=grd_mirror_r_template (grid, n, buffer, mask);
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_mirror_r (grid_t grid, int n, short *buffer, short mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=grd_mirror_r_template (grid, n, buffer, mask);
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_mirror (grid_t grid, int n, complex<float> *buffer, complex<float> mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=grd_mirror_r_template (grid, n, buffer, mask);
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_mirror_s (grid_t grid, int n, short *buffer, short mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float tmp;
  int i,j,m1,m2;

  for (j=0;j<grid.ny/2;j++) {
    for (i=0;i<grid.nx;i++) {
      m1=grid.nx*j+i;
      m2=grid.nx*(grid.ny-j-1)+i;
      tmp=buffer[m1];
      buffer[m1]=buffer[m2];
      buffer[m2]=tmp;
      }
    }
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int grd_extract_template (const char* filename, grid_t grid, int n, T *buf, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, ncid;
  int i1,i2,j,j1,j2;
  int ndimsp,var;
  size_t *start=NULL,*count=NULL;
  ptrdiff_t *stride=NULL;
  grid_t filegrid;

  status=grd_loadgrid(filename,&filegrid);

  status=nc_open(filename,0,&ncid);
  if (status !=NC_NOERR) {
    status=-1;
    return(status);
    }

/* grd is upside down */

/*   j1=NINT((grid.ymin-filegrid.ymin)/grid.dy); if grd was not upside down */
/*   j2=NINT((grid.ymax-filegrid.ymin)/grid.dy); */

  j1=int( NINT((filegrid.ymax-grid.ymax)/grid.dy) );
  j2=int( NINT((filegrid.ymax-grid.ymin)/grid.dy) );


  i1=int( NINT((grid.xmin-filegrid.xmin)/grid.dx) );
  i2=int( NINT((grid.xmax-filegrid.xmin)/grid.dx) );

  var=5;

  ndimsp=1;
  exitIfNull(
    start=(size_t *) malloc(ndimsp*sizeof(size_t))
    );
  exitIfNull(
    count=(size_t *) malloc(ndimsp*sizeof(size_t))
    );
  exitIfNull(
    stride=(ptrdiff_t *) malloc(ndimsp*sizeof(ptrdiff_t))
    );

  stride[0]=1;

  for (j=j1;j<=j2;j++) {
    start[0]=i1+j*filegrid.nx;
    count[0]=i2-i1+1;
    status=nc_get_vars_float(ncid,var,start,count,stride,&buf[(j-j1)*grid.nx]);
    if(status != NC_NOERR) goto error;
    }

  status=nc_close(ncid);
  
  free(start);
  free(count);
  free(stride);

  return(status);

  error:
  status=-1;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_extract (const char* filename, grid_t grid, int n, float *buf, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=grd_extract_template(filename, grid, n, buf, mask);
  return(status);
}
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_save(const char* filename, grid_t grid, int n, float *buf, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  historical float grd format, prior mirroring needed
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status, ncid;
  size_t lengthhp[2],index[10];
  int ndimsp,nvarsp,ngattsp,unlimdimidp,dim,var,att,length,dimids[50];
  char s1[1024]="side",s2[1024]="xysize";
  char *varname[6],*dimname[2], *attname[2], *attvalue[2];
  nc_type vartype[6];
  int varid[6],vardim[6],varndimsp[6];
  int i,idum;
  double ddum;
  float fdum;

  time_t creation_time;
  status=time(&creation_time);

  status=nc_create(filename,NC_NETCDF4,&ncid);
  if (status !=NC_NOERR) {
    status=-1;
    return(status);
    }

  ndimsp=2;

  dimname[0]=strdup(s1);
  lengthhp[0]=2;

  dimname[1]=strdup(s2);
  lengthhp[1]=(size_t) grid.nx*(size_t) grid.ny;

  nvarsp=6;
  ngattsp=2;
  unlimdimidp=-1;

  for (dim=0;dim<ndimsp;dim++) {
    status=nc_def_dim(ncid,dimname[dim],lengthhp[dim],&dimids[dim]);
    }

  varname[0]=(char *)("x_range");
  vartype[0]=NC_DOUBLE;
  varndimsp[0]= 1;
  vardim[0]= 0;
  varname[1]=(char *)("y_range");
  vartype[1]=NC_DOUBLE;
  varndimsp[1]= 1;
  vardim[1]= 0;
  varname[2]=(char *)("z_range");
  vartype[2]=NC_DOUBLE;
  varndimsp[2]= 1;
  vardim[2]= 0;
  varname[3]=(char *)("spacing");
  vartype[3]=NC_DOUBLE;
  varndimsp[3]= 1;
  vardim[3]= 0;
  varname[4]=(char *)("dimension");
  vartype[4]=NC_INT;
  varndimsp[4]= 1;
  vardim[4]= 0;
  varname[5]=(char *)("z");
  vartype[5]=NC_FLOAT;
  varndimsp[5]= 1;
  vardim[5]= 1;

  for (var=0;var<nvarsp;var++) {
    status=nc_def_var(ncid,varname[var],vartype[var],varndimsp[var],&vardim[var],&varid[var]);
    }

  ngattsp=2;
  attname[0]=strdup("title");
  attvalue[0]=strdup("GMT-like bathymetry file");
  att=0;
  length=strlen(attvalue[att])+1;
  length=1;
  status=nc_put_att_text(ncid,NC_GLOBAL,attname[att],length,attvalue[att]);

  att=1;
  attname[1]=strdup("source");
  exitIfNull(
    attvalue[1]=(char *) malloc(481)
    );
  for (i=0;i<481;i++) *(attvalue[1]+i)=0;
  strcpy(attvalue[1],"created by grd tools library, F. Lyard 2008");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,NC_GLOBAL,attname[att],length,attvalue[att]);

  attname[1]=strdup("creation");
  for (i=0;i<481;i++) *(attvalue[1]+i)=0;
  strcpy(attvalue[1],ctime(&creation_time));
  length=strlen(attvalue[att]);
  *(attvalue[1]+length-1)=0;
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,NC_GLOBAL,attname[att],length,attvalue[att]);

  att=0;
  attname[att]=strdup("units");
  attvalue[att]=strdup("user_x_unit");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,0,attname[att],length,attvalue[att]);

  free(attvalue[att]);
  attvalue[att]=strdup("user_y_unit");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,1,attname[att],length,attvalue[att]);

  free(attvalue[att]);
  attvalue[att]=strdup("user_z_unit");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,2,attname[att],length,attvalue[att]);

  free(attname[att]);
  free(attvalue[att]);
  attname[att]=strdup("long_name");
  attvalue[att]=(char *)("z");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,5,attname[att],length,attvalue[att]);

  free(attname[att]);
  attname[att]=strdup("scale_factor");
  fdum=1.;
  length=1;
  status=nc_put_att_float(ncid,5,attname[att],NC_FLOAT,length,&fdum);

  free(attname[att]);
  attname[att]=strdup("add_offset");
  fdum=0.;
  length=1;
  status=nc_put_att_float(ncid,5,attname[att],NC_FLOAT,length,&fdum);

  free(attname[att]);
  attname[att]=strdup("_FillValue");
  fdum=mask;
  length=1;
  status=nc_put_att_float(ncid,5,attname[att],NC_FLOAT,length,&fdum);

  free(attname[att]);
  attname[att]=strdup("node_offset");
  idum=0;
  length=1;
  status=nc_put_att_int(ncid,5,attname[att],NC_INT,length,&idum);

  status=nc_enddef(ncid);
  if(status != NC_NOERR) {
    nc_advise("nc_enddef", status, "ncid %d", ncid);
    return -1;
    }

//   status=nc_close(ncid);
//
//   status=nc_open(filename,NC_WRITE,&ncid);

  var=0;
  index[0]=0;
  ddum=grid.xmin;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  if(status != NC_NOERR) {
    nc_advise("ncattget", status, "ncid %d", ncid);
    return -1;
    }
  index[0]=1;
  ddum=grid.xmax;
  status=nc_put_var1_double(ncid,var,index,&ddum);

  var=1;
  index[0]=0;
  ddum=grid.ymin;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  index[0]=1;
  ddum=grid.ymax;
  status=nc_put_var1_double(ncid,var,index,&ddum);

  var=2;
  index[0]=0;
  ddum=0;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  index[0]=1;
  ddum=1;
  status=nc_put_var1_double(ncid,var,index,&ddum);

  var=3;
  index[0]=0;
  ddum=grid.dx;
  status=nc_put_var1_double(ncid,var,index,&ddum);
  index[0]=1;
  ddum=grid.dy;
  status=nc_put_var1_double(ncid,var,index,&ddum);

  var=4;
  index[0]=0;
  status=nc_put_var1_int(ncid,var,index,&(grid.nx));
  index[0]=1;
  status=nc_put_var1_int(ncid,var,index,&(grid.ny));

  var=5;
  status=nc_put_var_float(ncid,var,buf);

  status=nc_close(ncid);

  return(status);
  }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int grd_save(const char* filename, grid_t grid, int n, short *buf, short mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  historical short grd format, prior mirroring needed
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
  {
  int status, ncid;
  int var,att,length;
  char *varname[6],*attname[2], *attvalue[2];
  nc_type vartype[6];
  int varid[6],vardim[6],varndimsp[6];
  int idum;
  double ddum;
  short sdum;

  status=grd_savegrid(filename, grid);

  status=nc_open(filename,NC_WRITE,&ncid);
  if (status !=NC_NOERR) goto error;

  status=nc_redef(ncid);
  if (status !=NC_NOERR) goto error;

  var=5;
  varname[var]=(char *)("z");
  vartype[var]=NC_SHORT;
  varndimsp[var]= 1;
  vardim[var]= 1;

  status=nc_def_var(ncid,varname[var],vartype[var],varndimsp[var],&vardim[var],&varid[var]);
  if (status !=NC_NOERR) goto error;

  att=1;
  attname[att]=strdup("long_name");
  attvalue[att]=(char *)("\0");
  length=strlen(attvalue[att])+1;
  status=nc_put_att_text(ncid,5,attname[att],length,attvalue[att]);
  if(status != NC_NOERR) goto error;

  free(attname[att]);
  attname[att]=strdup("scale_factor");
  ddum=1.;
  length=1;
  status=nc_put_att_double(ncid,5,attname[att],NC_DOUBLE,length,&ddum);
  if(status != NC_NOERR) goto error;

  free(attname[att]);
  attname[att]=strdup("add_offset");
  ddum=0.;
  length=1;
  status=nc_put_att_double(ncid,5,attname[att],NC_DOUBLE,length,&ddum);
  if(status != NC_NOERR) goto error;

  free(attname[att]);
  attname[att]=strdup("_FillValue");
  sdum=mask;
  length=1;
  status=nc_put_att_short(ncid,5,attname[att],NC_SHORT,length,&sdum);

  free(attname[att]);
  attname[att]=strdup("node_offset");
  idum=0;
  length=1;
  status=nc_put_att_int(ncid,5,attname[att],NC_INT,length,&idum);
  if(status != NC_NOERR) goto error;

  status=nc_enddef(ncid);
  if (status !=NC_NOERR) goto error;

  var=5;
  status=nc_put_var_short(ncid,var,buf);
  if (status !=NC_NOERR) goto error;

  status=nc_close(ncid);
  if (status !=NC_NOERR) goto error;

  return(status);

  error:
  if(status != NC_NOERR) {
    nc_advise("ncattget", status, "ncid %d", ncid);
    return(status);
    }

  }

