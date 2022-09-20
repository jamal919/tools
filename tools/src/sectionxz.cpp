

/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Cyril Nguen        LA, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#define MAIN_SOURCE

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

 

#include "tools-structures.h"
     
#include "fe.h"
#include "map.h"
#include "netcdf-proto.h"
#include "archive.h"


   /* rank (number of dimensions) for each variable */
#  define RANK_X 2
#  define RANK_Z 2
#  define RANK_time 1
#  define RANK_T 3
#  define RANK_S 3
#  define RANK_U 3
#  define RANK_V 3
#  define RANK_rho 3

/*----------------------------------------------------------------------------*/
/**
\sa create_ncfile_xy() and create_ncfile_xyz()
*/
/*----------------------------------------------------------------------------*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int create_ncfile_xz(char *filename,size_t x_len, size_t z_len, grid_t grid,float mask)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/
{
   int  ncid; /* netCDF id */

   /* dimension ids */
   int x_dim;
   int z_dim;
   int time_counter_dim;

   /* dimension lengths */
/*
   size_t x_len = 100;
   size_t z_len = 100;
*/
   size_t time_counter_len = NC_UNLIMITED;

   /* variable ids */
   int X_id;
   int Z_id;
   int time_id;
   int T_id;
   int S_id;
   int U_id;
   int V_id;
   int rho_id;


   /* variable shapes */
   int X_dims[RANK_X];
   int Z_dims[RANK_Z];
   int time_dims[RANK_time];
   int T_dims[RANK_T];
   int S_dims[RANK_S];
   int U_dims[RANK_U];
   int V_dims[RANK_V];
   int rho_dims[RANK_rho];

   /* attribute vectors */
   float T_missing_value[1];
   float T__FillValue[1];
   float T_scale_factor[1];
   float T_add_offset[1];
   float S_missing_value[1];
   float S__FillValue[1];
   float S_scale_factor[1];
   float S_add_offset[1];
   float U_missing_value[1];
   float U__FillValue[1];
   float U_scale_factor[1];
   float U_add_offset[1];
   float V_missing_value[1];
   float V__FillValue[1];
   float V_scale_factor[1];
   float V_add_offset[1];
   float rho_missing_value[1];
   float rho__FillValue[1];
   float rho_scale_factor[1];
   float rho_add_offset[1];

   /* enter define mode */
   int stat = nc_create(filename, NC_CLOBBER, &ncid);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define dimensions */
   stat = nc_def_dim(ncid, "x", x_len, &x_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "z", z_len, &z_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "time_counter", time_counter_len, &time_counter_dim);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define variables */

   X_dims[0] = z_dim;
   X_dims[1] = x_dim;
   stat = nc_def_var(ncid, "X", NC_DOUBLE, RANK_X, X_dims, &X_id);
   nc_check_error(stat,__LINE__,__FILE__);

   Z_dims[0] = z_dim;
   Z_dims[1] = x_dim;
   stat = nc_def_var(ncid, "Z", NC_DOUBLE, RANK_Z, Z_dims, &Z_id);
   nc_check_error(stat,__LINE__,__FILE__);

   time_dims[0] = time_counter_dim;
   stat = nc_def_var(ncid, "time", NC_DOUBLE, RANK_time, time_dims, &time_id);
   nc_check_error(stat,__LINE__,__FILE__);

   T_dims[0] = time_counter_dim;
   T_dims[1] = z_dim;
   T_dims[2] = x_dim;
   stat = nc_def_var(ncid, "T", NC_FLOAT, RANK_T, T_dims, &T_id);
   nc_check_error(stat,__LINE__,__FILE__);

   S_dims[0] = time_counter_dim;
   S_dims[1] = z_dim;
   S_dims[2] = x_dim;
   stat = nc_def_var(ncid, "S", NC_FLOAT, RANK_S, S_dims, &S_id);
   nc_check_error(stat,__LINE__,__FILE__);

   U_dims[0] = time_counter_dim;
   U_dims[1] = z_dim;
   U_dims[2] = x_dim;
   stat = nc_def_var(ncid, "U", NC_FLOAT, RANK_U, U_dims, &U_id);
   nc_check_error(stat,__LINE__,__FILE__);

   V_dims[0] = time_counter_dim;
   V_dims[1] = z_dim;
   V_dims[2] = x_dim;
   stat = nc_def_var(ncid, "V", NC_FLOAT, RANK_V, V_dims, &V_id);
   nc_check_error(stat,__LINE__,__FILE__);

   rho_dims[0] = time_counter_dim;
   rho_dims[1] = z_dim;
   rho_dims[2] = x_dim;
   stat = nc_def_var(ncid, "rho", NC_FLOAT, RANK_rho, rho_dims, &rho_id);
   nc_check_error(stat,__LINE__,__FILE__);

   /* assign attributes */
   stat = nc_put_att_text(ncid, X_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, X_id, "long_name", 6, "X-axis");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, X_id, "standard_name", 6, "X-axis");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Z_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Z_id, "long_name", 5, "depth");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, Z_id, "standard_name", 5, "depth");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, time_id, "units", 7, "seconds");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, time_id, "calendar", 9, "gregorian");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, time_id, "title", 4, "Time");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, time_id, "long_name", 34, "Time axis elapsed from time_origin");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, time_id, "time_origin", 21, " 1950-JAN-01 00:00:00");
   nc_check_error(stat,__LINE__,__FILE__);
   T_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, T_id, "missing_value", NC_FLOAT, 1, T_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   T__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, T_id, "_FillValue", NC_FLOAT, 1, T__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   T_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, T_id, "scale_factor", NC_FLOAT, 1, T_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   T_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, T_id, "add_offset", NC_FLOAT, 1, T_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "standard_name", 21, "sea_water_temperature");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "long_name", 21, "sea water temperature");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "units", 2, "\260C");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "content", 3, "TZX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, T_id, "associate", 12, "time z x");
   nc_check_error(stat,__LINE__,__FILE__);
   S_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, S_id, "missing_value", NC_FLOAT, 1, S_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   S__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, S_id, "_FillValue", NC_FLOAT, 1, S__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   S_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, S_id, "scale_factor", NC_FLOAT, 1, S_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   S_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, S_id, "add_offset", NC_FLOAT, 1, S_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "standard_name", 18, "sea_water_salinity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "long_name", 18, "sea water salinity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "units", 6, "PSU   ");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "content", 3, "TZX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, S_id, "associate", 12, "time z x");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "units", 5, "m s-1");
   nc_check_error(stat,__LINE__,__FILE__);
   U_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, U_id, "missing_value", NC_FLOAT, 1, U_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   U__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, U_id, "_FillValue", NC_FLOAT, 1, U__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   U_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, U_id, "scale_factor", NC_FLOAT, 1, U_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   U_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, U_id, "add_offset", NC_FLOAT, 1, U_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "long_name", 28, "ocean barotropic E-component");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "standard_name", 27, "eastward_sea_water_velocity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "content", 3, "TZX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, U_id, "associate", 12, "time z x");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "units", 5, "m s-1");
   nc_check_error(stat,__LINE__,__FILE__);
   V_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, V_id, "missing_value", NC_FLOAT, 1, V_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   V__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, V_id, "_FillValue", NC_FLOAT, 1, V__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   V_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, V_id, "scale_factor", NC_FLOAT, 1, V_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   V_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, V_id, "add_offset", NC_FLOAT, 1, V_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "long_name", 22, "ocean barotropic E-component");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "standard_name", 25, "northward_sea_water_velocity");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "content", 3, "TZX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, V_id, "associate", 12, "time z x");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "units", 3, "PSI");
   nc_check_error(stat,__LINE__,__FILE__);
   rho_missing_value[0] = mask;
   stat = nc_put_att_float(ncid, rho_id, "missing_value", NC_FLOAT, 1, rho_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   rho__FillValue[0] = mask;
   stat = nc_put_att_float(ncid, rho_id, "_FillValue", NC_FLOAT, 1, rho__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   rho_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, rho_id, "scale_factor", NC_FLOAT, 1, rho_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   rho_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, rho_id, "add_offset", NC_FLOAT, 1, rho_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "long_name", 13, "water density");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "standard_name", 25, "water_density");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "content", 3, "TZX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, rho_id, "associate", 12, "time z x");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 5, "CF1.0");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", 7, "regular");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", 7, "whatsoever.nc");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "production", 32, "Florent Lyard produced this file");
   nc_check_error(stat,__LINE__,__FILE__);

   /* leave define mode */
   stat = nc_enddef (ncid);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_close(ncid);
   nc_check_error(stat,__LINE__,__FILE__);
   return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float *rbuf;
  float *rbufx,*rbufy;
  double t0,dT,pulsation;
  double *serie[500],a1,p1,a2,p2,zr,zi,d;

  float x,y;
  float  dummy,rmask,*rbuffer;
  float  *rbufferx,*rbuffery;
  float  spec=-9999;
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int nitems;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL,*sectionfile=NULL,*rootname=NULL;
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*path=NULL;
  grid_t grid;
  mesh_t mesh;
  int nwave=0,column,ncid,iterative=0,count;
 
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

        case 's' :
          sscanf(argv[n+1],"%f",&spec);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          rootname= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
          sectionfile= strdup(argv[n]);
          n++;
          break;
        break;
      }
      free(keyword);
    }

  rmask=spec;

 if(path ==NULL) path=strdup(".");

 if(sectionfile==NULL) {
   iterative=1;
   sectionfile=(char *)malloc(1024);
   sprintf(sectionfile,"%s001.xz",rootname);
   }

 if(meshfile != NULL) {
    bool debug=false;
    status=fe_readmesh_QUODDY(meshfile, sectionfile, 0, &mesh, debug);
    if(status !=0) goto error;
/*
    status=fe_list(&mesh);
    if(status !=0) goto error;
*/
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }
  mesh.type=1;
  nndes=mesh.nvtxs;
  //fe_initaffine(&mesh);

  grid=get_grid_n(mesh);
  status=map_completegridaxis_2(&grid);

  elts=fe_scan_elements(mesh,grid,1);
  if(elts==NULL) goto error;

  output=(char *)malloc(1024);

  sprintf(output,"%s.xz.nc",rootname);

  rmask=-9999.;
  status=create_ncfile_xz(output, grid.nx, grid.ny, grid,rmask);
  if(status !=0) goto error;

  status=nc_open(output,NC_WRITE,&ncid);
  if(status !=0) goto error;

/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=nc_writeaxis(ncid,grid);
  if(status !=0) goto error;
   
  rbuffer=(float *)malloc(nndes*sizeof(float));
  rbuf=(float *)malloc(grid.nx*grid.ny*sizeof(fcomplex));

  count=0;
  do
  {
  count++;
  sprintf(sectionfile,"%s%3.3d.xz",rootname,count);
  printf("reading %s\n",sectionfile);
  for (column=0;column<5;column++) {
    status=section_loadr1(sectionfile, mesh.nvtxs, rbuffer,column);
    if(status!=0) {
      printf("cannot read %s\n",sectionfile);
      goto finished;
      }
    status=fe_map(mesh,rbuffer,grid,elts,rbuf,rmask);
    if(status!=0) goto error;
/**----------------------------------------------------------------------------
    obsolete call, will be suppressed */
    status=nc_write_r1(ncid, grid, column+3, count-1, rbuf,rmask);
    if(status!=0) goto error;
    }
  if(iterative==0) break;
  } while(status==0);
  
 finished:
  free(rbuffer);
  free(rbuf);

  status = nc_close(ncid);

  printf("end of gridit ... \n");
  free(elts);
  __ERR_BASE_LINE__("exiting\n");exit(0);
error:
 __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
