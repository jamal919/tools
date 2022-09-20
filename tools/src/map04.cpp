
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

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <stdarg.h>
#include <map>

#include "tools-structures.h"

#include "map.def"
#include "map.h"
#include "grd.h"
#include "ascii.h"
#include "netcdf-proto.h"
#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "map.def"



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_loadfield_grd_template(const char *input, grid_t *grid, T **buffer, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,m,n,status;

  status=grd_loadgrid(input,grid);
  if(status !=0) {
    __OUT_BASE_LINE__("cannot load bathymetry file=%s\n",input);
    exit(-1);
    }

  *buffer= new T[grid->nx*grid->ny];
  status= grd_loadr1(input,*grid,*buffer,mask);

  if(isnan(*mask)) {
    *mask=1.e+10;
    *mask=99999;
/* *----------------------------------------------------------------------------
    nan is very inconvenient as mask value*/
    for(m=0;m<grid->nx*grid->ny;m++) {
      if(isnan((*buffer)[m])) {
        (*buffer)[m]=*mask;
        }
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield_grd(const char *input, grid_t *grid, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_grd_template(input, grid, buffer, mask);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield_grd(const char *input, grid_t *grid, short **buffer, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_grd_template(input, grid, buffer, mask);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_loadfield_raster_template(const char *input, grid_t *grid, T*  & buffer, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status;
  char tile='g';
  int    nitems,count=0;
  FILE   *in;
  short *buf2;


  grid->dx=1./120.;
  grid->dy=1./120.;
  
  grid->nx=21600/2;
  
  tile=input[0];
  
  switch (tile) {
    case 'a':
      grid->xmin=-180;
      grid->ymin=50;
      grid->ny=4800;
      break;
    case 'b':
      grid->xmin=-90;
      grid->ymin=50;
      grid->ny=4800;
      break;
    case 'c':
      grid->xmin=0;
      grid->ymin=50;
      grid->ny=4800;
      break;
    case 'd':
      grid->xmin=90;
      grid->ymin=50;
      grid->ny=4800;
      break;
    case 'e':
      grid->xmin=-180;
      grid->ymin=0;
      grid->ny=6000;
      break;
    case 'f':
      grid->xmin=-90;
      grid->ymin=0;
      grid->ny=6000;
      break;
    case 'g':
      grid->xmin=0;
      grid->ymin=0;
      grid->ny=6000;
      break;
    case 'h':
      grid->xmin=90;
      grid->ymin=0;
      grid->ny=6000;
      break;
    case 'i':
      grid->xmin=-180;
      grid->ymin=-50;
      grid->ny=6000;
      break;
    case 'j':
      grid->xmin=-90;
      grid->ymin=-50;
      grid->ny=6000;
      break;
    case 'k':
      grid->xmin=0;
      grid->ymin=-50;
      grid->ny=6000;
      break;
    case 'l':
      grid->xmin=90;
      grid->ymin=-50;
      grid->ny=6000;
      break;
    case 'm':
      grid->xmin=-180;
      grid->ymin=-90;
      grid->ny=4800;
      break;
    case 'n':
      grid->xmin=-90;
      grid->ymin=-90;
      grid->ny=4800;
      break;
    case 'o':
      grid->xmin=0;
      grid->ymin=-90;
      grid->ny=4800;
      break;
    case 'p':
      grid->xmin=90;
      grid->ymin=-90;
      grid->ny=4800;
      break;
    }
  
  grid->xmax=grid->xmin+(grid->nx-1)*grid->dx;
  grid->ymax=grid->ymin+(grid->ny-1)*grid->dy;
  
  
  grid->modeH=0;
  
  status=map_completegridaxis(grid,2);
  
  status=-1;

  buf2=new short[grid->nx*grid->ny];
  buffer=new T[grid->nx*grid->ny];

  in=fopen(input,"r");

read:
  nitems=fread(buf2,2,grid->nx*grid->ny,in);
  if(nitems!=grid->nx*grid->ny) goto error;

  for (j=0;j<grid->ny;j++) {
    for (i=0;i<grid->nx;i++) {
      k=j*grid->nx+i;
      l=(grid->ny-j-1)*grid->nx+i;
      buffer[l]=buf2[k];
      if(buffer[k]==0.) buffer[k]=1e+35;
      }
    }

  fclose(in);
  status=0;
 error:
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield_raster(const char *input, grid_t *grid, float*  & buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_raster_template(input, grid, buffer, mask);
  return(status);
}

  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield_raster(const char *input, grid_t *grid, short*  & buffer, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_raster_template(input, grid, buffer, mask);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_loadfield_cdf_template(const char *input,const char *varname, const char *xname, const char *yname, grid_t *grid, T **buffer, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,n,v,status;
  int    verbose,frame=0;
  cdfgbl_t data_info,grid_info;
  cdfvar_t var_info, x_info, y_info;

  verbose=0;
  status=cdf_globalinfo(input,&data_info,verbose);
  for (v=0;v<data_info.nvarsp;v++) {
    printf("variable %3d: name %s, type %d,ndim %d \n",v,(data_info.variable[v]).name,data_info.variable[v].type,data_info.variable[v].ndim);
    }
  status=cdf_varinfo(input,varname,&var_info,verbose);
  if(status!=0) return(-1);

  status= cdf_globalinfo(input,&grid_info,verbose);
  status= poc_getgrid2d (input, grid_info, var_info, grid);

  *buffer=new T[grid->nx*grid->ny];
  status= poc_getvar2d (input, var_info.id, frame, *buffer, mask, var_info);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_loadfield_cdf_new_template(const char *input,const char *varname, const char *xname, const char *yname, grid_t *grid, T **buffer, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,n,v,status;
  int    verbose,frame=0;
  cdfgbl_t data_info,grid_info;
  cdfvar_t var_info, x_info, y_info;
  poc_var_t var;
  T default_mask=1.e+36;
  double scale, offset;

  verbose=0;
  
  status=poc_get_grid (input, varname, grid);
  if(status!=0) return(status);
  
  status=poc_inq_var(input, varname, &var,0);
  if(status!=0) return(status);

  *buffer=new T[grid->nx*grid->ny];
  status=poc_get_var(input, varname, *buffer, verbose);
  if(status!=0) return(status);
  
  status=poc_decode_mask(var,&scale,&offset,mask,default_mask);

  status=0;
  
  poc_scale_data(*buffer,grid->Hsize(),scale,offset,mask,verbose);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield_cdf(const char *input,const char *varname, const char *xname, const char *yname, grid_t *grid, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_cdf_new_template(input, varname, xname, yname, grid, buffer, mask);
  return(status);
}

  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield_cdf(const char *input,const char *varname, const char *xname, const char *yname, grid_t *grid, short **buffer, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_cdf_template(input, varname, xname, yname, grid, buffer, mask);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield_cdf(const char *input, const char *varname, grid_t & grid, float* & buffer, float & mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=map_loadfield_cdf_template(input, varname, 0, 0, &grid, &buffer, &mask);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int map_loadfield_template(const char *filename, int format, const char *varname, const char *xname, const char *yname, const char *proj, grid_t *grid, T **buffer, T *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,m,n,status=-1;
  string s="";

  if(filename==0) {
    check_error(-1, "input file name not defined", __LINE__, __FILE__, 1);
    }
  s.assign(filename);
  
  switch(format) {
    case MAP_FORMAT_CODE_GRD:
      status=map_loadfield_grd(filename,grid, buffer, mask);
      break;
    case MAP_FORMAT_CODE_COMODO:
      if(varname==0) {
        printf("netcdf varname not specified, abort\n");
        return(-1);
        }
      else {
        status=map_loadfield_cdf(filename, varname, xname, yname, grid, buffer, mask);
        }
      break;
    case MAP_FORMAT_CODE_GLOBE:
      status=map_loadfield_raster(filename, grid, *buffer, mask);
      break;
    case MAP_FORMAT_CODE_ASCII:
      status=ascii_load(filename, grid, buffer, mask);
      break;
    case MAP_FORMAT_CODE_ASCII_GIS:
      status=ascii_load_GIS(filename, grid, proj, buffer, mask);
      proj=0;
      break;
    case MAP_FORMAT_CODE_ASCII_SLIM:
      status=ascii_load_SLIM(filename, grid, proj, buffer, mask);
      proj=0;
      break;
    default:
      if(varname==0) {
        printf("netcdf varname not specified, abort\n");
        return(-1);
        }
      status=map_loadfield_cdf(filename, varname, xname, yname, grid, buffer, mask);
      break;
    }
    
  if(proj!=0) {
    if(grid->modeH==0) status=map_completegridaxis(grid,2);
    status=map_projection_backward(*grid, proj);
    }

    
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield(const char *filename, int format, const char *varname, const char *xname, const char *yname, const char *proj, grid_t *grid, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_template(filename, format, varname, xname, yname, proj, grid, buffer, mask);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield(const char *filename, int format, const char *varname, grid_t *grid, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_template(filename, format, varname, 0, 0, 0, grid, buffer, mask);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield(const char *filename, int format, const char *varname, const char *xname, const char *yname, const char *proj, grid_t *grid, short **buffer, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_template(filename, format, varname, xname, yname, proj, grid, buffer, mask);
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield(const char *filename, int format, const char *varname, grid_t *grid, short **buffer, short *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=map_loadfield_template(filename, format, varname, 0, 0, 0, grid, buffer, mask);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_extractfield(const char *filename,grid_t grid, float **buffer, float *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,m,n,status;
  string s;

  s.assign(filename);

  *buffer=new float[grid.nx*grid.ny];

  if(s.find(".grd")!=-1) {
    status=grd_extract (filename, grid, grid.nx, *buffer, mask);
    }
  else if(s.find(".nc")!=-1) {
//    status=map_loadfield_cdf(filename,grid, buffer, mask);
    }
  return(0);
}

