
/***************************************************************************
T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative
  
E-mail: florent.lyard@legos.obs-mip.fr
***************************************************************************/
/**
\file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Yves Soufflet      LEGOS, Toulouse, France
**/


#include "version-macros.def" //for VERSION and REVISION

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <algorithm>

#include <errno.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "fe.h"
#include "poc-time.h"
#include "list.h"
#include "map.h"
#include "map.def"
#include "topo.h"
#include "grd.h"
#include "geo.h"
#include "cefmo.h"
#include "functions.h"
#include "tides.h"
#include "netcdf-proto.h"
#include "poc-netcdf-data.hpp"
#include "sturm-liouville.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield3D(const char *datafile, const char *gridfile, const char *var, int frame, grid_t & grid, float * &buffer, float & mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,n,v,status;
  int    vdepth;
  cdfgbl_t g_info,grid_info;
  cdfvar_t v_info;
  
  buffer=0;

  status=cdf_globalinfo(datafile,&g_info,verbose);
  if(status!=0) return(-1);
  
  if(verbose) {
    for (v=0;v<g_info.nvarsp;v++) {
      printf("variable %3d: name %s, type %d,ndim %d \n",v,(g_info.variable[v]).name,g_info.variable[v].type,g_info.variable[v].ndim);
      }
    }
    
  status=cdf_varinfo(datafile,var,&v_info,verbose);
  if(status!=0) return(-1);

  status= poc_getgrid3d (gridfile, g_info, v_info, &grid, &vdepth);
  if(status!=0) return(-1);

  buffer=new float[grid.nx*grid.ny*grid.nz];
  status= poc_getvar3d (datafile, v_info.id, frame, buffer, &mask, v_info);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_loadfield3D(const char *datafile, const char *gridfile,const char *var, grid_t & grid, float * &buffer, float & mask, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    status;
  int    frame=0;

  status=map_loadfield3D(datafile, gridfile, var, frame, grid, buffer, mask,verbose);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int save_SGXY_C_template(const char *output, grid_t & grid, complex<T> *buffer, complex<T> mask, 
                                        const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  fcomplex z;
  cdfvar_t variable, variableA, variableG;
  pocgrd_t ncgrid;
  date_t origin;
  
  printf("#################################################################\n");
  printf("write %s %s in output file : %s\n",varnameA,varnameG, output);
  if(createfile==1) {
    status= poc_createfile(output);
    }
    
  if(creategrid==1) {
    status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
    }
  else {
    status=poc_inq_grid(output,"","",&ncgrid);
    }

  T spec=(T) abs(mask);
  
  poc_standardvariable_xy(&variableA, varnameA,spec,unit,(T) 1., (T) 0.,standard,standard,standard,ncgrid);
  status=create_ncvariable(output, &variableA);
  
  poc_standardvariable_xy(&variableG,varnameG,spec,"degrees",(T) 1., (T) 0.,standard,standard,standard,ncgrid);
  status=create_ncvariable(output, &variableG);
  
  T* tmp=new T[grid.Hsize()];
  
  for(size_t n=0;n<grid.Hsize();n++) {
    if(buffer[n]!=mask) {
      tmp[n]=abs(buffer[n]);
      if(!isnormal(tmp[n]) and tmp[n]!=0) {
        printf("%s anomaly\n",__FUNCTION__);
        tmp[n]=spec;
        }
      }
    else {
      tmp[n]=spec;
      }
    }
 status=poc_write_xy(output, grid, variableA.id, tmp);
 
  for(size_t n=0;n<grid.Hsize();n++) {
    if(buffer[n]!=mask) {
      tmp[n]=-arg(buffer[n])*180./M_PI;
      if(!isnormal(tmp[n]) and tmp[n]!=0) {
        printf("%s anomaly\n",__FUNCTION__);
        tmp[n]=spec;
        }
      }
    else {
      tmp[n]=spec;
      }
    }
  status=poc_write_xy(output, grid, variableG.id, tmp);
  
  delete[] tmp;
  
  variableA.destroy();
  variableG.destroy();
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXY_C(const char *output, grid_t & grid, complex<float> *buffer, complex<float> mask, 
                                        const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  
  status=save_SGXY_C_template(output, grid, buffer, mask, varnameA, varnameG, unit, standard, createfile, creategrid);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXY_C(const char *output, grid_t & grid, complex<double> *buffer, complex<double> mask, 
                                        const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  
  status=save_SGXY_C_template(output, grid, buffer, mask, varnameA, varnameG, unit, standard, createfile, creategrid);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int save_SGXYZ_template(const char *output, grid_t & grid, T *buffer, T mask, const char *varname,const char *unit, const char *standard, int createfile, int creategrid, const char *hlabel, const char *vlabel)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;

  
  printf("#################################################################\n");
  printf("write %s in output file : %s\n",varname, output);
  if(createfile==1) {
    status= poc_createfile(output);
    }

  if(creategrid==1) {
    status=poc_sphericalgrid_xyz(output,hlabel,vlabel,grid,&ncgrid);
    }
  else {
    status=poc_inq_grid(output,hlabel,vlabel,&ncgrid);
    }

  status=poc_standardvariable_xyz(varname,mask,unit,1., 0.,standard,standard,standard,ncgrid, variable);
  status=create_ncvariable(output, &variable);
  status=poc_write_xyz(output,  grid, variable.id,buffer);
  variable.destroy();

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYZ(const char *output, grid_t grid, float *buffer, float mask, const char *varname,const char *unit, const char *standard, int createfile, int creategrid, const char *hlabel, const char *vlabel)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,status;

  status=save_SGXYZ_template(output, grid, buffer, mask, varname, unit,  standard, createfile, creategrid, hlabel, vlabel);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYZ(const char *output, grid_t grid, double *buffer, double mask, const char *varname,const char *unit, const char *standard, int createfile, int creategrid, const char *hlabel, const char *vlabel)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,status;

  status=save_SGXYZ_template(output, grid, buffer, mask, varname, unit,  standard, createfile, creategrid, hlabel, vlabel);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int save_SGXYZ_C_template(const char *output, grid_t & grid, complex<T> *buffer, complex<T> mask, 
                                        const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid, const char *hlabel, const char *vlabel)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  fcomplex z;
  cdfvar_t variable, variableA, variableG;
  pocgrd_t ncgrid;
  date_t origin;
  
  printf("#################################################################\n");
  printf("write %s %s in output file : %s\n",varnameA,varnameG, output);
  if(createfile==1) {
    status= poc_createfile(output);
    }
    
  if(creategrid==1) {
    status=poc_sphericalgrid_xyz(output,hlabel,vlabel,grid,&ncgrid);
    }
  else {
    status=poc_inq_grid(output,hlabel,vlabel,&ncgrid);
    }

  T spec=(T) abs(mask);
  
  status=poc_standardvariable_xyz(varnameA,spec,unit,(T) 1., (T) 0.,standard,standard,standard,ncgrid,variableA);
  status=create_ncvariable(output, &variableA);
  
  status=poc_standardvariable_xyz(varnameG,spec,"degrees",(T) 1., (T) 0.,standard,standard,standard,ncgrid,variableG);
  status=create_ncvariable(output, &variableG);
  
  T* tmp=new T[grid.size()];
  
  for(size_t n=0;n<grid.size();n++) {
    if(buffer[n]!=mask) {
      tmp[n]=abs(buffer[n]);
      if(!isnormal(tmp[n]) and tmp[n]!=0) {
        printf("%s anomaly\n",__FUNCTION__);
        tmp[n]=spec;
        }
      }
    else {
      tmp[n]=spec;
      }
    }
 status=poc_write_xyz(output, grid, variableA.id, tmp);
  for(size_t n=0;n<grid.size();n++) {
    if(buffer[n]!=mask) {
      tmp[n]=-arg(buffer[n])*180./M_PI;
      if(!isnormal(tmp[n]) and tmp[n]!=0) {
        printf("%s anomaly\n",__FUNCTION__);
        tmp[n]=spec;
        }
      }
    else {
      tmp[n]=spec;
      }
    }
  status=poc_write_xyz(output, grid, variableG.id, tmp);
  
  delete[] tmp;
  
  variableA.destroy();
  variableG.destroy();
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYZ_C(const char *output, grid_t & grid, complex<float> *buffer, complex<float> mask, 
                                        const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid, const char *hlabel, const char *vlabel)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  
  status=save_SGXYZ_C_template(output, grid, buffer, mask, varnameA, varnameG, unit, standard, createfile, creategrid, hlabel, vlabel);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYZ_C(const char *output, grid_t & grid, complex<double> *buffer, complex<double> mask, 
                                        const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid, const char *hlabel, const char *vlabel)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  
  status=save_SGXYZ_C_template(output, grid, buffer, mask, varnameA, varnameG, unit, standard, createfile, creategrid, hlabel, vlabel);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int save_SGXYT_template(const char *output, grid_t grid, T *buffer, T mask, const char *varname,const char *unit, const char *standard, 
                                       int createfile, int creategrid, int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  date_t origin;
  
  printf("#################################################################\n");
  printf("write %s (frame %d) in output file : %s\n",varname, frame, output);
  if(createfile==1) {
    status= poc_createfile(output);
    }
    
  if(creategrid==1) {
    grid.nz=1;
    grid.z=NULL;
    status=poc_sphericalgrid_xyt(output,"",grid,&ncgrid);
    if(grid.time!=0) {
      variable= poc_standardtime("time","nt","seconds", origin);
      status=create_ncvariable(output, &variable);
      for(count=0; count<grid.nt; count++) status=poc_writetime(output, count, variable.id,  grid.time[count]);
      variable.destroy();
      }
    }
  else {
    status=poc_inq_grid(output,"","",&ncgrid);
    }

  status=poc_standardvariable_xyt(varname,mask,unit,(T)1.,(T)0.,standard,standard,standard,ncgrid,variable);
  status=create_ncvariable(output, &variable);
  
  status=poc_write_xyt(output,  grid, frame, variable.id, buffer);

  variable.destroy();
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYT(const char *output, grid_t & grid, float *buffer, float mask, const char *varname,const char *unit, const char *standard,  int createfile, int creategrid, int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  
  status=save_SGXYT_template(output, grid, buffer, mask,  varname, unit, standard, createfile, creategrid, frame);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYT(const char *output, grid_t & grid, double *buffer, double mask, const char *varname,const char *unit, const char *standard,  int createfile, int creategrid, int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  
  status=save_SGXYT_template(output, grid, buffer, mask,  varname, unit, standard, createfile, creategrid, frame);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int save_SGXYT_template(const char *output, grid_t grid, T **buffer, T mask, const char *varname,const char *unit, const char *standard, 
                                       int createfile, int creategrid, int range[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  date_t origin;
  int backup=grid.nz;
  
  printf("#################################################################\n");
  printf("write %s (%d to %d) in output file : %s\n",varname, range[0], range[1], output);
  if(createfile==1) {
    status= poc_createfile(output);
    }
    
  if(creategrid==1) {
    grid.nz=1;
    grid.z=NULL;
    status=poc_sphericalgrid_xyt(output,"",grid,&ncgrid);
    if(grid.time!=0) {
      variable= poc_standardtime("time","nt","seconds", origin);
      status=create_ncvariable(output, &variable);
      for(count=0; count<grid.nt; count++) status=poc_writetime(output, count, variable.id,  grid.time[count]);
      variable.destroy();
      }
    }
  else {
    status=poc_inq_grid(output,"","",&ncgrid);
    }

  status=poc_standardvariable_xyt(varname,mask,unit,(T)1.,(T)0.,standard,standard,standard,ncgrid,variable);
  status=create_ncvariable(output, &variable);

  size_t frame=0;
  for(frame=range[0];frame<=range[1];frame++) {
    status=poc_write_xyt(output,  grid, frame, variable.id, buffer[frame]);
    }
  variable.destroy();
  grid.nz=backup;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYT(const char *output, grid_t & grid, float **buffer, float mask, const char *varname,const char *unit, const char *standard,
		 int createfile, int creategrid, int range[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  
  status=save_SGXYT_template(output, grid, buffer, mask,  varname, unit, standard, createfile, creategrid, range);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYT(const char *output, grid_t & grid, double **buffer, double mask, const char *varname,const char *unit, const char *standard,
		 int createfile, int creategrid, int range[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  
  status=save_SGXYT_template(output, grid, buffer, mask,  varname, unit, standard, createfile, creategrid, range);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYT_C(const char *output, grid_t grid, complex<double> **buffer, complex<double> mask, 
                   const char *varnameA, const char *varnameG, const char *unit, const char *standard, int createfile, int creategrid, int range[2])
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  fcomplex z;
  cdfvar_t variable, variableA, variableG;
  pocgrd_t ncgrid;
  date_t origin;
  
  printf("#################################################################\n");
  printf("write %s %s (%d to %d) in output file : %s\n",varnameA,varnameG, range[0], range[1], output);
  if(createfile==1) {
    status= poc_createfile(output);
    }
    
  if(creategrid==1) {
    grid.nz=1;
    grid.z=NULL;
    status=poc_sphericalgrid_xyt(output,"",grid,&ncgrid);
    if(grid.time!=0) {
      variable= poc_standardtime("time","t","seconds", origin);
      status=create_ncvariable(output, &variable);
      for(count=0; count<grid.nt; count++) status=poc_writetime(output, count, variable.id,  grid.time[count]);
      variable.destroy();
      }
    }
  else {
//     status=poc_inq_grid(output,hlabel,vlabel,&ncgrid);
    }

  double spec=(double) abs(mask);
  
  status=poc_standardvariable_xyt(varnameA,spec,unit,(double) 1., (double) 0.,standard,standard,standard,ncgrid,variableA);
  status=create_ncvariable(output, &variableA);
  
  status=poc_standardvariable_xyt(varnameG,spec,"degrees",(double) 1., (double) 0.,standard,standard,standard,ncgrid,variableG);
  status=create_ncvariable(output, &variableG);
  
  double* tmp=new double[grid.Hsize()];
  
  size_t frame=0;
  for(frame=range[0];frame<=range[1];frame++) {
    for(size_t n=0;n<grid.Hsize();n++) {
      if(buffer[frame][n]!=mask) {
        tmp[n]=abs(buffer[frame][n]);
        if(!isnormal(tmp[n]) and tmp[n]!=0) {
          printf("%s anomaly\n",__FUNCTION__);
          tmp[n]=spec;
          }
        }
      else {
        tmp[n]=spec;
        }
      }
    range_t<double> r=poc_minmax(tmp,grid.Hsize(),spec);
    status=poc_write_xyt(output,  grid, frame, variableA.id, tmp);
    for(size_t n=0;n<grid.Hsize();n++) {
      if(buffer[frame][n]!=mask) {
        tmp[n]=-arg(buffer[frame][n])*180./M_PI;;
        if(!isnormal(tmp[n]) and tmp[n]!=0) {
          printf("%s anomaly\n",__FUNCTION__);
          tmp[n]=spec;
          }
        }
      else {
        tmp[n]=spec;
        }
      }
    status=poc_write_xyt(output,  grid, frame, variableG.id, tmp);
    }
  
  delete[] tmp;
  
  variableA.destroy();
  variableG.destroy();
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYZT_C(const char *output, grid_t grid, complex<double> *buffer, complex<double> mask, const char *varnameA, const char *varnameG, const char *unit, const char *standard, 
                   int createfile, int creategrid, const char * vlabel, int frame)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,status;
  fcomplex z;
  cdfvar_t variable, variableA, variableG;
  pocgrd_t ncgrid;
  date_t origin;
  
  
  printf("#################################################################\n");
  printf("write %s %s (frame %d) in output file : %s\n",varnameA,varnameG, frame, output);
  if(createfile==1) {
    status= poc_createfile(output);
    }
/*------------------------------------------------------------------------
  check wether file exist or not*/
  FILE *in = fopen(output, "r");
  if(in == NULL) {
    status= poc_createfile(output);
    }
  else {
    fclose(in);
    }

    
  if(creategrid==1) {
    status=poc_sphericalgrid_xyzt(output,"",vlabel,grid,&ncgrid);
    if(grid.time!=0) {
      variable= poc_standardtime("time","t","seconds", origin);
      status=create_ncvariable(output, &variable);
      for(count=0; count<grid.nt; count++) status=poc_writetime(output, count, variable.id,  grid.time[count]);
      variable.destroy();
      }
    }
  else {
    status=poc_inq_grid(output,"",vlabel,&ncgrid);
    }

  double spec=(double) abs(mask);
  
  status=poc_standardvariable_xyzt(varnameA,spec,unit,(double) 1., (double) 0.,standard,standard,standard,ncgrid,variableA);
  status=create_ncvariable(output, &variableA);
  
  status=poc_standardvariable_xyzt(varnameG,spec,"degrees",(double) 1., (double) 0.,standard,standard,standard,ncgrid,variableG);
  status=create_ncvariable(output, &variableG);
  
  double* tmp=new double[grid.size()];
  
  for(size_t n=0;n<grid.size();n++) {
    if(buffer[n]!=mask) {
      tmp[n]=abs(buffer[n]);
      if(!isnormal(tmp[n]) and tmp[n]!=0) {
        printf("%s anomaly\n",__FUNCTION__);
        tmp[n]=spec;
        }
      }
    else {
      tmp[n]=spec;
      }
    }
    
  status=poc_write_xyzt(output,  grid, frame, variableA.id, tmp);
  
  for(size_t n=0;n<grid.size();n++) {
    if(buffer[n]!=mask) {
      tmp[n]=-arg(buffer[n])*180./M_PI;;
      if(!isnormal(tmp[n]) and tmp[n]!=0) {
        printf("%s anomaly\n",__FUNCTION__);
        tmp[n]=spec;
        }
      }
    else {
      tmp[n]=spec;
      }
    }
  status=poc_write_xyzt(output,  grid, frame, variableG.id, tmp);
  
  delete[] tmp;
  
  variableA.destroy();
  variableG.destroy();
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYZT(const char *output, grid_t grid, float **buffer, float mask, const char *varname,const char *unit, const char *standard, 
		  int createfile, int creategrid, const char *hlabel, const char *vlabel, int range[2])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  float *rbufx=NULL,*rbufy=NULL,*topo=NULL, rmask;

  
  printf("#################################################################\n");
  printf("write %s (%d to %d) in output file : %s\n",varname, range[0], range[1], output);
  if(createfile==1) {
    status= poc_createfile(output);
    }
  if(creategrid==1) {
    status=poc_sphericalgrid_xyzt(output,hlabel,vlabel,grid,&ncgrid);
    }
  else {
    status=poc_inq_grid(output,hlabel,vlabel,&ncgrid);
    }

  status=poc_standardvariable_xyzt(varname,mask,unit,1., 0.,standard,standard,standard,ncgrid, variable);
  
  status=create_ncvariable(output, &variable);
  
  size_t frame=0;
  for(frame=range[0];frame<=range[1];frame++) {
    status=poc_write_xyzt(output,  grid, frame, variable.id, buffer[frame]);
    }
  variable.destroy();

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int save_SGXYZT_template(const char *output, grid_t grid, T *buffer, T mask, const char *varname, const char *unit, const char *standard, 
                        int createfile, int creategrid, const char *hlabel, const char *vlabel, int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,status;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  
  printf("#################################################################\n");
  printf("write %s, frame=%d in output file : %s\n",varname, frame, output);
  if(createfile==1) {
    status= poc_createfile(output);
    }
    
  if(creategrid==1) {
    status=poc_sphericalgrid_xyzt(output,hlabel,vlabel,grid,&ncgrid);
    }
  else {
    status=poc_inq_grid(output,hlabel,vlabel,&ncgrid);
    }

  status=poc_standardvariable_xyzt(varname,mask,unit,1., 0.,standard,standard,standard,ncgrid, variable);
  
  status=create_ncvariable(output, &variable);
  
  status=poc_write_xyzt(output,  grid, frame, variable.id, buffer);

  variable.destroy();
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_SGXYZT(const char *output, grid_t grid, float *buffer, float mask, const char *varname, const char *unit, const char *standard, 
                        int createfile, int creategrid, const char *hlabel, const char *vlabel, int frame)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,status;
  
  status=save_SGXYZT_template(output,grid, buffer, mask, varname, unit, standard, createfile, creategrid, hlabel, vlabel, frame);
  
  return(status);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int VerticalModes_SpectralDecomposition_backup(int maxmodes, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*------------------------------------------------------------------------------
// 
//   Interface to decompose harmonic fields into modal contributions
//   
//   Presently hard-coded for NEMO-derived files, badly-formed structured 
//   ouputs are just a pain in the ... as usual
//   
//     v f v f v f   
//     t u t u t u   last T line unusable as v (above) masked  
//     v f v f v f
//     t u t u t u   first T line unusable as v start after
//     
//   
//   u,v and t have same dimensions: 2 lines and columns too much in t
// 
// ------------------------------------------------------------------------------*/
// {
//   int i,k,l,n,status;
//   char *modesfile=0,*ufile=0;
//   grid_t grid,ugrid,vgrid;
//   double **modes=0, **Umodes=0, **Pmodes=0, mask;
//   int MinUmodes;
//   complex<double> **u=0, **v=0, umask, *buffer=0, *u_buffer=0, *v_buffer=0, *ssh=0, cmask;
//   complex<double> **u_decomposition=0, **v_decomposition=0, **p_decomposition=0;
//   complex<double> *ubarotropic=0, *ubaroclinic=0, *vbarotropic=0, *vbaroclinic=0, *pbarotropic=0, *pbaroclinic=0;
//   string filename,varnameA,varnameG,varname;
//   int nlayers;
//   int verbose=0;
//   poc_data_t<float> scalarData;
//  
// /*------------------------------------------------------------------------------
//   tidal currents */
// 
// /*------------------------------------------------------------------------------
//   load velocity vertical modes */
//   filename="NEMO.vertical-modes.nc";
//   varname="Umodes";
//   
//   status=poc_get_grid(filename, varname, &grid, verbose, -1);
//   
//   MinUmodes=grid.nt;
//   nlayers=grid.nz;
//   
//   MinUmodes=nlayers;
//   
//   int nmasked;
//   status=scalarData.init(filename,varname);
//   for(i=0;i<MinUmodes;i++) {  
//     status=scalarData.read_data(filename,i);
//     MinUmodes=min(scalarData.nframes, maxmodes);
//     if(Umodes==0) Umodes=new double*[MinUmodes];
//     Umodes[i]=new double[grid.size()];
//     for(n=0;n<grid.size();n++) Umodes[i][n]=scalarData.data[n];
//     mask=scalarData.mask;
//     nmasked=occurence<double>(mask, Umodes[i], grid.size());
//     }
//     
//   umask=NC_FILL_COMPLEX;
//   u=new complex<double>* [grid.Hsize()];
//   v=new complex<double>* [grid.Hsize()];
//   for(n=0;n<grid.Hsize();n++) {
//     u[n]=new complex<double>[nlayers];
//     v[n]=new complex<double>[nlayers];
//     for(k=0;k<nlayers;k++) {
//       u[n][k]=umask;
//       v[n][k]=umask;
//       }
//     }
//   
// /*------------------------------------------------------------------------------
//   zonal velocity, load and interpolate at tracer nodes */
//   filename="M2-vozocrtx-atlas.nc";
//   varnameA="vozocrtx_a";
//   varnameG="vozocrtx_G";
//   
//   status=poc_get_grid(filename, varnameA, &ugrid, verbose, -1);
//   
//   buffer=new complex<double>[ugrid.size()];
//   
//   status=poc_get_cvara(filename,varnameA, varnameG,-1, buffer);
//   
//   poc_data_t<int8_t> t_mask,u_mask,v_mask;
//   poc_data_t<float> e3t;
//   status=t_mask.read_data("mesh_mask.nc","tmask",-1,1);
//   status=u_mask.read_data("mesh_mask.nc","umask",-1,1);
//   status=v_mask.read_data("mesh_mask.nc","vmask",-1,1);
//   status=e3t.read_data("mesh_mask.nc","e3t",-1,1);
//     
//   for(int j=2;j<grid.ny-2;j++) {
//     for(i=2;i<grid.nx-2;i++) {
//       size_t n=j*grid.nx+i;
//       for(k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+j*grid.nx+i;
//         if(t_mask.data[m]==0) continue;
//         if(Umodes[0][m]==mask) {
//           continue;
//           }
//         size_t m1=k*ugrid.Hsize()+j*ugrid.nx+i-1;
//         size_t m2=k*ugrid.Hsize()+j*ugrid.nx+i;
//         if(u_mask.data[m1]==0 and t_mask.data[m]==1) buffer[m1]=0;
//         if(u_mask.data[m2]==0 and t_mask.data[m]==1) buffer[m2]=0;
//         if(buffer[m1]!=umask and buffer[m2]!=umask) {
//           u[n][k]=0.5*(buffer[m1]+buffer[m2]);
//           }
//         else  {
// //           u[n][k]=umask;
//           u[n][k]=0.0;
//           }
//         }
//       }
//     }
//   
//   delete[] buffer;
//   ugrid.free();
//   
// /*------------------------------------------------------------------------------
//   meridian velocity, load and interpolate at tracer nodes */
//   filename="M2-vomecrty-atlas.nc";
//   varnameA="vomecrty_a";
//   varnameG="vomecrty_G";
//   
//   status=poc_get_grid(filename, varnameA, &vgrid, verbose, -1);
//   
//   buffer=new complex<double>[vgrid.size()];
//   
//   status=poc_get_cvara(filename,varnameA, varnameG,-1, buffer);
//   
//   for(int j=2;j<grid.ny-2;j++) {
//     for(int i=2;i<grid.nx-2;i++) {
//       size_t n=j*grid.nx+i;
//       for(k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+j*grid.nx+i;
//         if(t_mask.data[m]==0) continue;
//         if(Umodes[0][m]==mask) {
//           continue;
//           }
//         size_t m1=k*vgrid.Hsize()+(j-1)*vgrid.nx+i;
//         size_t m2=k*vgrid.Hsize()+j*vgrid.nx+i;
//         if(v_mask.data[m1]==0 and t_mask.data[m]==1) buffer[m1]=0;
//         if(v_mask.data[m2]==0 and t_mask.data[m]==1) buffer[m2]=0;
//         if(buffer[m1]!=umask and buffer[m2]!=umask) {
//           v[n][k]=0.5*(buffer[m1]+buffer[m2]);
//           }
//         else  {
//           v[n][k]=0;
// //           v[n][k]=umask;
//           }
//         }
//       }
//     }
//   delete[] buffer;
//   vgrid.free();
//     
//   u_decomposition=new complex<double>*[MinUmodes];
//   for(i=0;i<MinUmodes;i++) {
//     u_decomposition[i]=new complex<double>[grid.Hsize()];
//     for(n=0;n<grid.Hsize();n++) u_decomposition[i][n]=umask;
//     }
//   v_decomposition=new complex<double>*[MinUmodes];
//   for(i=0;i<MinUmodes;i++) {
//     v_decomposition[i]=new complex<double>[grid.Hsize()];
//     for(n=0;n<grid.Hsize();n++) v_decomposition[i][n]=umask;
//     }
//     
//   ubarotropic=new complex<double>[grid.size()];
//   for(size_t m=0;m<grid.size();m++) ubarotropic[m]=umask;
//   ubaroclinic=new complex<double>[grid.size()];
//   for(size_t m=0;m<grid.size();m++) ubaroclinic[m]=umask;
//   
//   vbarotropic=new complex<double>[grid.size()];
//   for(size_t m=0;m<grid.size();m++) vbarotropic[m]=umask;
//   vbaroclinic=new complex<double>[grid.size()];
//   for(size_t m=0;m<grid.size();m++) vbaroclinic[m]=umask;
//     
//   int nRequestedProcs=-1;
//   int nprocs=initialize_OPENMP(nRequestedProcs);
//   
// #pragma omp parallel for private(status) if(nprocs>1)
//   for(int jj=2;jj<grid.ny-2;jj++) {
//     complex<double> *decomposition=new complex<double>[nlayers]; 
//     double **lmodes=new double*[MinUmodes];
//     for(int i=0;i<MinUmodes;i++) {
//       lmodes[i]=new double[nlayers];
//       }
//     for(int ii=2;ii<grid.nx-2;ii++) {
//       size_t n=jj*grid.nx+ii;
//       int nvalidlayers=0;
//       int ntruemodes=0;
//       int count;
//       bool skip=false;
//       for(int k=0;k<nlayers;k++) {
//         if(u[n][k]!=umask) nvalidlayers++;
//         }
//       for(int i=0;i<MinUmodes;i++) {
//         double std=0;
//         count=0;
//         for(int k=0;k<nlayers;k++) {
//           size_t m=k*grid.Hsize()+n;
//           lmodes[i][k]=Umodes[i][m];
//           if(lmodes[i][k]!=mask) {
//             std+=lmodes[i][k]*lmodes[i][k];
//             count++;
//             }
//           }
//         if(count!=0 and std!=0.0) {
// /*------------------------------------------------------------------------------
//           count effective number of valid */        
// 	  ntruemodes++;
// /*------------------------------------------------------------------------------
//           check added because of tracer and velocity grids mismatches */        
//           if(count!=nvalidlayers) {
//             skip=true;
//             }
// /*------------------------------------------------------------------------------
//           normalize vertical modes */        
//           std=sqrt(std/(double) count);
//           for(int k=0;k<nlayers;k++) {
//             size_t m=k*grid.Hsize()+n;
//             if(lmodes[i][k]!=mask) lmodes[i][k]/=std;
//             if(Umodes[i][m]!=mask) Umodes[i][m]/=std;
//             }
//           }
//         else {
//           if(i==0) skip=true;
//           }
//         }
// /*------------------------------------------------------------------------------
//       check added because of tracer and velocity grids mismatches */
//       int nmodes=min(ntruemodes,MinUmodes);
//       if(skip) continue;
//       status=Umodes_decomposition_v2(nvalidlayers, lmodes, nmodes, u[n], decomposition, debug, 0);
//       for(int i=0;i<MinUmodes;i++) {
//         u_decomposition[i][n]=decomposition[i];
//         }
//       for(int k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+n;
//         if(Umodes[0][m]!=mask) ubarotropic[m]=lmodes[0][k]*decomposition[0];
//         ubaroclinic[m]=0;
//         for(int i=1;i<nmodes;i++) {
//           if(Umodes[i][m]!=mask) ubaroclinic[m]+=lmodes[i][k]*decomposition[i];
//           double chk=abs(ubaroclinic[m]);
//           if(!isnormal(chk) and chk!=0) {
//             printf("%s anomaly\n",__FUNCTION__);
//             }
//           }
//         }
//       status=Umodes_decomposition_v2(nvalidlayers, lmodes, nmodes, v[n], decomposition, debug, 0);
//       for(int i=0;i<MinUmodes;i++) {
//         v_decomposition[i][n]=decomposition[i];
//         }
//       for(int k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+n;
//         if(Umodes[0][m]!=mask) vbarotropic[m]=lmodes[0][k]*decomposition[0];
//         vbaroclinic[m]=0;
//         for(int i=1;i<nmodes;i++) {
//           if(Umodes[i][m]!=mask) vbaroclinic[m]+=lmodes[i][k]*decomposition[i];
//           double chk=abs(ubaroclinic[m]);
//           if(!isnormal(chk) and chk!=0) {
//             printf("%s anomaly\n",__FUNCTION__);
//             }
//           }
//         }
//       }
//     delete[] decomposition;
//     deletep2D(&lmodes,MinUmodes);
//     }
//     
//   int range[2]={0,MinUmodes-1};
//   
//   int create_file=1;
//   int create_grid=1;
//   
// /*------------------------------------------------------------------------------
//   modal decomposition coefficients */
//   status=save_SGXYT_C("modal-decomposition.nc", grid, u_decomposition, umask, "Uc_a", "Uc_G", "dimensionless", "u_coeffcients", create_file, create_grid, range);
//   create_file=0;
//   create_grid=1;
//   status=save_SGXYT_C("modal-decomposition.nc", grid, v_decomposition, umask, "Vc_a", "Vc_G", "dimensionless", "u_coeffcients", create_file, create_grid, range);
// 
//   complex<double> *laplacian=0;
//   int ref=0;
//   float *wavelength=new float[grid.Hsize()], rmask=1.e+10;
//     
// //   for(int i=0;i<MinUmodes;i++) {
//   status=map_laplacian(grid, u_decomposition[1], umask, ref, laplacian);
//   for(size_t m=0;m<grid.Hsize();m++) {
//     wavelength[m]=rmask;
//     if(laplacian[m]!=umask and u_decomposition[1][m]!=0.0) {
//       complex<double> ratio=laplacian[m]/u_decomposition[1][m];
//       double k2=-real(ratio);
//       if(k2>0.0) wavelength[m]=2.*M_PI/sqrt(k2);
//       else {
//         printf("%s anomaly\n",__FUNCTION__);
//         }
//       }
//     }
//   create_file=1;
//   status=save_SGXY("modal-diagnostics.nc", grid, wavelength, rmask, "uL", "m", "uL", create_file);
// //     }
//   
//   create_file=0;
//   create_grid=1;
// /*------------------------------------------------------------------------------
//   modal fields */
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, ubarotropic, umask, "Ubt_a", "Ubt_G", "m/s", "u_barotropic", create_file, create_grid);
//   create_file=0;
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, ubaroclinic, umask, "Ubc_a", "Ubc_G", "m/s", "u_baroclinic", create_file, create_grid);
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, vbarotropic, umask, "Vbt_a", "Vbt_G", "m/s", "v_barotropic", create_file, create_grid);
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, vbaroclinic, umask, "Vbc_a", "Vbc_G", "m/s", "v_baroclinic", create_file, create_grid);
//     
//   int nkept=5;
//   for(i=0;i<nkept;i++) {
//     for(n=0;n<grid.Hsize();n++) {
//       for(int k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+n;
//         ubaroclinic[m]=umask;
//         vbaroclinic[m]=umask;
//         if(Umodes[i][m]==mask) continue;
//         ubaroclinic[m]=Umodes[i][m]*u_decomposition[i][n];
//         vbaroclinic[m]=Umodes[i][m]*v_decomposition[i][n];
//         }
//       }
//     create_file=(i==0);
//     create_grid=(i==0);
//     int frame=i;
//     status=save_SGXYZT_C("recomposed.nc", grid, ubaroclinic, umask, "u_modal_a","u_modal_G","m/s", "u_modal", create_file, create_grid,frame);
//     create_file=0;
//     create_grid=0;
//     status=save_SGXYZT_C("recomposed.nc", grid, vbaroclinic, umask, "v_modal_a","v_modal_G","m/s", "v_modal", create_file, create_grid,frame);
//     }
//       
// /*------------------------------------------------------------------------------
//   pressure */
//   filename="NEMO.vertical-modes.nc";
//   varname="Pmodes";
//   
//   status=scalarData.init(filename,varname);
//   for(i=0;i<MinUmodes;i++) {  
//     status=scalarData.read_data(filename,i);
//     MinUmodes=min(scalarData.nframes, maxmodes);
//     if(Pmodes==0) {
//       Pmodes=new double*[MinUmodes];
//       for(k=0;k<MinUmodes;k++) Pmodes[k]=0;
//       }
//     if(Pmodes[i]==0) Pmodes[i]=new double[grid.size()];
//     for(n=0;n<grid.size();n++) Pmodes[i][n]=scalarData.data[n];
//     mask=scalarData.mask;
//     nmasked=occurence<double>(mask, Pmodes[i], grid.size());
//     }
//     
//   filename="M2-Pbc.nc";
//   varnameA="Pbc_a";
//   varnameG="Pbc_G";
// 
//   buffer=new complex<double>[grid.size()];
//   
//   status=poc_get_cvara(filename,varnameA, varnameG,-1, buffer);
//   
//   umask=NC_FILL_COMPLEX;
//   
//   for(n=0;n<grid.Hsize();n++) {
//     for(k=0;k<nlayers;k++) {
//       size_t m=k*grid.Hsize()+n;
//       if(t_mask.data[m]==0) continue;
//       if(Pmodes[0][m]==mask) {
//         continue;
//         }
//       if(buffer[m]!=umask) {
//         u[n][k]=buffer[m];
//         }
//       else  {
//         u[n][k]=umask;
//         }
//       }
//     }
// 
// #if 0
//   filename="M2-sossheig-atlas.nc";
//   varnameA="sossheig_a";
//   varnameG="sossheig_G";
//   
//   status=poc_get_cvara(filename,varnameA, varnameG,-1, buffer);
//   
//   umask=NC_FILL_COMPLEX;
// //   u=new complex<double>* [grid.Hsize()];
//   
//   for(n=0;n<grid.Hsize();n++) {
// //     u[n]=new complex<double>[nlayers];
//     if(buffer[n]!=umask) {
//       buffer[n]*=1019.2*9.81;
//       for(k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+n;
//         if(u[n][k]!=umask) {
//           u[n][k]+=buffer[n];
//           }
//         }
//       }
//     }
// #endif
// 
//   delete[] buffer;
//   
//   p_decomposition=new complex<double>*[MinUmodes];
//   for(i=0;i<MinUmodes;i++) {
//     p_decomposition[i]=new complex<double>[grid.Hsize()];
//     for(n=0;n<grid.Hsize();n++) p_decomposition[i][n]=umask;
//     }
//     
//   pbarotropic=new complex<double>[grid.size()];
//   for(size_t m=0;m<grid.size();m++) pbarotropic[m]=umask;
//   pbaroclinic=new complex<double>[grid.size()];
//   for(size_t m=0;m<grid.size();m++) pbaroclinic[m]=umask;
//     
// //   for(n=0;n<grid.Hsize();n++) {
// #pragma omp parallel for private(status) if(nprocs>1)
//   for(int jj=2;jj<grid.ny-2;jj++) {
//     complex<double> *decomposition=new complex<double>[nlayers]; 
//     double **lmodes=new double*[MinUmodes];
//     for(int i=0;i<MinUmodes;i++) {
//       lmodes[i]=new double[nlayers];
//       }
//     for(int ii=2;ii<grid.nx-2;ii++) {
//       size_t n=jj*grid.nx+ii;
//     int ntruelayers=0;
//     int ntruemodes=0;
//     double count;
//     bool skip=false;
//     for(int k=0;k<nlayers;k++) {
//       if(u[n][k]!=umask) ntruelayers++;
//       }
//     for(int i=0;i<MinUmodes;i++) {
//       double std=0;
//       count=0;
//       for(int k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+n;
//         lmodes[i][k]=Pmodes[i][m];
//         if(lmodes[i][k]!=mask) {
//           std+=lmodes[i][k]*lmodes[i][k];
//           count+=1.0;
//           }
//         }
//       if(count!=0 and std!=0.0) {
// /*------------------------------------------------------------------------------
//         count effective number of valid */        
// 	ntruemodes++;
// /*------------------------------------------------------------------------------
//         check added because of tracer and velocity grids mismatches */        
//         if(count!=ntruelayers) {
//           skip=true;
//           }
//         std=sqrt(std/count);
//         for(int k=0;k<nlayers;k++) {
//           size_t m=k*grid.Hsize()+n;
//           if(lmodes[i][k]!=mask) lmodes[i][k]/=std;
//           if(Pmodes[i][m]!=mask) Pmodes[i][m]/=std;
//           }
//         }
//       else {
//         if(i==0) skip=true;
//         }
//       }
// //     debug=true;
// /*------------------------------------------------------------------------------
//     check added because of tracer and velocity grids mismatches */
//     int nmodes;
//     if(MinUmodes>ntruelayers) {
// //       skip=true;
//       nmodes=min(ntruemodes,ntruelayers);
//       }
//     else {
//       nmodes=min(ntruemodes,MinUmodes);
//       }
//     if(skip) continue;
//     status=Umodes_decomposition_v2(ntruelayers, lmodes, nmodes, u[n], decomposition, debug, 0);
//     for(int i=0;i<MinUmodes;i++) {
//       p_decomposition[i][n]=decomposition[i];
//       }
//     for(int k=0;k<nlayers;k++) {
//       size_t m=k*grid.Hsize()+n;
//       if(Pmodes[0][m]!=mask) pbarotropic[m]=lmodes[0][k]*decomposition[0];
//       pbaroclinic[m]=0;
//       for(int i=1;i<nmodes;i++) {
//         if(Pmodes[i][m]!=mask) pbaroclinic[m]+=lmodes[i][k]*decomposition[i];
//         double chk=abs(pbaroclinic[m]);
//         if(!isnormal(chk) and chk!=0) {
//           printf("%s anomaly\n",__FUNCTION__);
//           }
//         }
//       }
//     }
//     deletep2D(&lmodes,MinUmodes);
//     delete[] decomposition;
//     }
//     
//   create_file=0;
//   create_grid=1;
// /*------------------------------------------------------------------------------
//   modal decomposition coefficients */
//   status=save_SGXYT_C("modal-decomposition.nc", grid, p_decomposition, umask, "Pc_a", "Pc_G", "dimensionless", "u_coeffcients", create_file, create_grid, range);
// 
//   status=map_laplacian(grid, p_decomposition[1], umask, ref, laplacian);
//   for(size_t m=0;m<grid.Hsize();m++) {
//     wavelength[m]=rmask;
//     if(laplacian[m]!=umask and p_decomposition[1][m]!=0.0) {
//       complex<double> ratio=laplacian[m]/p_decomposition[1][m];
//       double k2=-real(ratio);
//       if(k2>0.0) wavelength[m]=2.*M_PI/sqrt(k2);
//       else {
//         printf("%s anomaly\n",__FUNCTION__);
//         }
//       }
//     }
//   create_file=0;
//   status=save_SGXY("modal-diagnostics.nc", grid, wavelength, rmask, "pL", "m", "pL", create_file);
//   
//   create_file=0;
//   create_grid=1;
// /*------------------------------------------------------------------------------
//   modal fields */
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, pbarotropic, umask, "Pbt_a", "Pbt_G", "N/m²", "p_barotropic", create_file, create_grid);
//   create_file=0;
//   create_grid=1;
//   status=save_SGXYZ_C("modal-decomposition.nc", grid, pbaroclinic, umask, "Pbc_a", "Pbc_G", "N/m²", "p_baroclinic", create_file, create_grid);
//   
//   for(i=0;i<nkept;i++) {
//     for(n=0;n<grid.Hsize();n++) {
//       for(int k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+n;
//         pbaroclinic[m]=umask;
//         if(Pmodes[i][m]==mask) continue;
//         pbaroclinic[m]=Pmodes[i][m]*p_decomposition[i][n];
//         }
//       }
//     create_file=0;
//     create_grid=0;
//     int frame=i;
//     status=save_SGXYZT_C("recomposed.nc", grid, pbaroclinic, umask, "p_modal_a","p_modal_G","N/m²", "p_modal", create_file, create_grid,frame);
//     }
// 
//   
//   double *uflux,*vflux,dmask=NC_FILL_DOUBLE;
//   uflux=new double[grid.Hsize()];
//   vflux=new double[grid.Hsize()];
//   
//   for(n=0;n<grid.Hsize();n++) {
//     uflux[n]=0;
//     vflux[n]=0;
//     double weight=0;
//     for(k=0;k<nlayers;k++) {
//       size_t m=k*grid.Hsize()+n;
//       if(u[n][k]!=umask) {
//         uflux[n]+=e3t.data[m]*integrale_time_2(pbarotropic[m], ubarotropic[m]);
//         vflux[n]+=e3t.data[m]*integrale_time_2(pbarotropic[m], vbarotropic[m]);
//         weight+=e3t.data[m];
//         }
//       }
//     if(weight==0) {
//       uflux[n]=dmask;
//       vflux[n]=dmask;
//       }
//     }
//   create_file=1;
//   create_grid=1;
// /*------------------------------------------------------------------------------
//   modal fields */
//   status=save_SGXY("modal-energy.nc", grid, uflux, dmask, "uFbt_x", "N/s", "barotropic_energy_flux", 1);
//   status=save_SGXY("modal-energy.nc", grid, vflux, dmask, "uFbt_y", "N/s", "barotropic_energy_flux", 0);
//   
//   for(n=0;n<grid.Hsize();n++) {
//     uflux[n]=0;
//     vflux[n]=0;
//     double weight=0;
//     for(k=0;k<nlayers;k++) {
//       size_t m=k*grid.Hsize()+n;
//       if(u[n][k]!=umask) {
// 	uflux[n]+=e3t.data[m]*integrale_time_2(pbaroclinic[m], ubaroclinic[m]);
//         vflux[n]+=e3t.data[m]*integrale_time_2(pbaroclinic[m], vbaroclinic[m]);
// 	weight+=e3t.data[m];
//         }
//       }
//     if(weight==0) {
//       uflux[n]=dmask;
//       vflux[n]=dmask;
//       }
//     }
//   create_file=0;
//   create_grid=0;
// /*------------------------------------------------------------------------------
//   modal fields */
//   status=save_SGXY("modal-energy.nc", grid, uflux, dmask, "uFbc_x", "N/s", "baroclinic_energy_flux", 0);
//   status=save_SGXY("modal-energy.nc", grid, vflux, dmask, "uFbc_y", "N/s", "baroclinic_energy_flux", 0);
//   
//   for(i=0;i<nkept;i++) {
//     for(n=0;n<grid.Hsize();n++) {
//       for(int k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+n;
//         pbaroclinic[m]=umask;
//         ubaroclinic[m]=umask;
//         vbaroclinic[m]=umask;
//         if(Pmodes[i][m]==mask) continue;
//         pbaroclinic[m]=Pmodes[i][m]*p_decomposition[i][n];
//         ubaroclinic[m]=Umodes[i][m]*u_decomposition[i][n];
//         vbaroclinic[m]=Umodes[i][m]*v_decomposition[i][n];
//         }
//       uflux[n]=0;
//       vflux[n]=0;
//       double weight=0;
//       for(k=0;k<nlayers;k++) {
//         size_t m=k*grid.Hsize()+n;
//         if(pbaroclinic[m]!=umask) {
// 	  uflux[n]+=e3t.data[m]*integrale_time_2(pbaroclinic[m], ubaroclinic[m]);
//           vflux[n]+=e3t.data[m]*integrale_time_2(pbaroclinic[m], vbaroclinic[m]);
// 	  weight+=e3t.data[m];
//           }
//         }
//       if(weight==0) {
//         uflux[n]=dmask;
//         vflux[n]=dmask;
//         }
//       }
//     create_file=0;
//     create_grid=(i==0);
//     int frame=i;
//     status=save_SGXYT("modal-energy.nc", grid, uflux, dmask, "Fmodal_x","N/m²", "Fmodal_x", create_file, create_grid, frame);
//     create_grid=0;
//     status=save_SGXYT("modal-energy.nc", grid, vflux, dmask, "Fmodal_y","N/m²", "Fmodal_y", create_file, create_grid, frame);
//     }
//     
//   deletep2D(&p_decomposition,MinUmodes);
//   deletep2D(&u_decomposition,MinUmodes);
//   deletep2D(&v_decomposition,MinUmodes);
// 
// /*------------------------------------------------------------------------------
//   only for verification, commented to keep files small */ 
// //   for(n=0;n<grid.Hsize();n++) {
// //     for(k=0;k<nlayers;k++) {
// //       size_t m=k*grid.Hsize()+n;
// //       if(u[n][k]!=umask) pbarotropic[m]=u[n][k]/9.81/1019.2;
// //       else pbarotropic[m]=umask;
// //       }
// //     }
// //   create_file=0;
// //   create_grid=1;
// //   status=save_SGXYZ_C("modal-decomposition.nc", grid, pbarotropic, umask, "P_a", "P_G", "m", "u_coeffcients", create_file, create_grid);
//  
//   return(0);    
// }
