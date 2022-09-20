
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

/*----------------------------------------------------------------------------*/
/**
\sa create_ncfile_xz() and create_ncfile_xyz()
*/
/*----------------------------------------------------------------------------*/
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

int create_ncfile_xy(char *filename,size_t x_len, size_t y_len, grid_t grid,float mask)

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/
{		

   int  ncid;			/* netCDF id */

   /* dimension ids */
   int x_dim;
   int y_dim;
   int time_counter_dim;

   /* dimension lengths */
/*
   size_t x_len = 100;
   size_t y_len = 100;
*/
   size_t time_counter_len = NC_UNLIMITED;

   /* variable ids */
   int lon_id;
   int lat_id;
   int time_id;
   int h_id;

   /* rank (number of dimensions) for each variable */
#  define RANK_lon 2
#  define RANK_lat 2
#  define RANK_time 1
#  define RANK_h 3

   /* variable shapes */
   int lon_dims[RANK_lon];
   int lat_dims[RANK_lat];
   int time_dims[RANK_time];
   int h_dims[RANK_h];

   /* attribute vectors */
   double lon_valid_min[1];
   double lon_valid_max[1];
   double lat_valid_min[1];
   double lat_valid_max[1];
   float h_missing_value[1];
   float h__FillValue[1];
   float h_scale_factor[1];
   float h_add_offset[1];

   /* enter define mode */
   int stat = nc_create(filename, NC_CLOBBER, &ncid);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define dimensions */
   stat = nc_def_dim(ncid, "x", x_len, &x_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "y", y_len, &y_dim);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_def_dim(ncid, "time_counter", time_counter_len, &time_counter_dim);
   nc_check_error(stat,__LINE__,__FILE__);

   /* define variables */

   lon_dims[0] = y_dim;
   lon_dims[1] = x_dim;
   stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
   nc_check_error(stat,__LINE__,__FILE__);

   lat_dims[0] = y_dim;
   lat_dims[1] = x_dim;
   stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
   nc_check_error(stat,__LINE__,__FILE__);

   time_dims[0] = time_counter_dim;
   stat = nc_def_var(ncid, "time", NC_DOUBLE, RANK_time, time_dims, &time_id);
   nc_check_error(stat,__LINE__,__FILE__);

   h_dims[0] = time_counter_dim;
   h_dims[1] = y_dim;
   h_dims[2] = x_dim;
   stat = nc_def_var(ncid, "h", NC_FLOAT, RANK_h, h_dims, &h_id);
   nc_check_error(stat,__LINE__,__FILE__);

   /* assign attributes */
   stat = nc_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
   nc_check_error(stat,__LINE__,__FILE__);
   lon_valid_min[0] = -180;
   stat = nc_put_att_double(ncid, lon_id, "valid_min", NC_DOUBLE, 1, lon_valid_min);
   nc_check_error(stat,__LINE__,__FILE__);
   lon_valid_max[0] = 180;
   stat = nc_put_att_double(ncid, lon_id, "valid_max", NC_DOUBLE, 1, lon_valid_max);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lon_id, "nav_model", 12, "Default grid");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
   nc_check_error(stat,__LINE__,__FILE__);
   lat_valid_min[0] = -90;
   stat = nc_put_att_double(ncid, lat_id, "valid_min", NC_DOUBLE, 1, lat_valid_min);
   nc_check_error(stat,__LINE__,__FILE__);
   lat_valid_max[0] = 90;
   stat = nc_put_att_double(ncid, lat_id, "valid_max", NC_DOUBLE, 1, lat_valid_max);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, lat_id, "nav_model", 12, "Default grid");
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
   stat = nc_put_att_text(ncid, h_id, "units", 1, "m");
   nc_check_error(stat,__LINE__,__FILE__);
   h_missing_value[0] = -9999.;
   stat = nc_put_att_float(ncid, h_id, "missing_value", NC_FLOAT, 1, h_missing_value);
   nc_check_error(stat,__LINE__,__FILE__);
   h__FillValue[0] = -9999.;
   stat = nc_put_att_float(ncid, h_id, "_FillValue", NC_FLOAT, 1, h__FillValue);
   nc_check_error(stat,__LINE__,__FILE__);
   h_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, h_id, "scale_factor", NC_FLOAT, 1, h_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__);
   h_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, h_id, "add_offset", NC_FLOAT, 1, h_add_offset);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, h_id, "long_name", 19, "sea level elevation");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, h_id, "standard_name", 21, "sea_surface_elevation");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, h_id, "float_name", 1, "h");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, h_id, "content", 3, "TYX");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, h_id, "associate", 12, "time lat lon");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 5, "CF1.0");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", 7, "regular");
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", 13, "whatsoever.nc");
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

  grid_t get_zonegrid_fatal(char *zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  grid_t grid;
  int zone_initialised=0;

/* NEA grid */
  if(strcmp(zone,"NEA")==0) {
  map_set2Dgrid(&grid,-20.0,+29.5,+15.0,+65.0,1./10.,1./10.);
  zone_initialised=1;
  }

/* NEA_hr grid haute resolution */
  if(strcmp(zone,"NEA_hr")==0) {
  map_set2Dgrid(&grid,-20.0,+29.5,+15.0,+65.0,1./50.,1./50.);
  zone_initialised=1;
  }

/* Manche_hr grid haute resolution */
  if(strcmp(zone,"Manche_hr")==0) {
  map_set2Dgrid(&grid,-6.0,+6.0,+48.0,+53.0,1./100.,1./100.);
  zone_initialised=1;
  }

 /* NEA_shom grid distribution shom  haute resolution */
  if(strcmp(zone,"NEA_shom")==0) {
  map_set2Dgrid(&grid,-20.0,+40.0,+15.0,+55.0,1./30.,1./30.);
  zone_initialised=1;
  }
      
/* Iroise grid */
  if(strcmp(zone,"iroise")==0) {
  map_set2Dgrid(&grid,-6.0,+47.25,-3.5,+49.0,1./60.,1./60.);
  zone_initialised=1;
  }
    
/* global grid */
  if(strcmp(zone,"global")==0) {
  map_set2Dgrid(&grid,0.0,-80.0,+360.0,+80.0,1./4.,1./4.);
  zone_initialised=1;
  }
  
/* global grid */
  if(strcmp(zone,"global-loren")==0) {
  map_set2Dgrid(&grid,0.0,-90.0,+360.0,+90.0,1.,1.);
  zone_initialised=1;
  }
  
    
/* mesdsea grid, normal resolution */
  if(strcmp(zone,"medsea")==0) {
  map_set2Dgrid(&grid,-10.0,27.5,40.0,47.5,1./10.,1./10.);
  zone_initialised=1;
  }

/* alobran grid, normal resolution */
  if(strcmp(zone,"alboran")==0) {
  map_set2Dgrid(&grid,-9.0,32,-2+1./300,38+1./300,1./60.,1./60.);
  zone_initialised=1;
  }
  
/* alobran grid, normal resolution */
  if(strcmp(zone,"strait-of-sicily")==0) {
  map_set2Dgrid(&grid,-9.0,32,-2+1./300,38+1./300,1./60.,1./60.);
  zone_initialised=1;
  }
 
  if(zone_initialised!=1) {
    __OUT_BASE_LINE__("no valid region specified (%s); abort...\n",zone);
    exit(-1);
    }

  return(grid);

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
  char *meshfile=NULL,*depthfile=NULL,*output=NULL,*path=NULL,*nodefile=NULL;
  grid_t grid;
  mesh_t mesh;
  int nwave=0,column,ncid,iterative=0,count;
 
  fct_echo(argc,argv);
 
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          nodefile= strdup(argv[n+1]);
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

        case 'z' :
          zone= strdup(argv[n+1]);
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
   }
 
 if(meshfile != NULL) {
    status=fe_readmesh_TGL(meshfile,nodefile,&mesh);
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

  if(zone !=NULL) grid=get_zonegrid_fatal(zone);
  else            grid=get_grid_n(mesh);
 
  status=map_completegridaxis_2(&grid);

  elts=fe_scan_elements(mesh,grid,mesh.type);
  if(elts==NULL) goto error;

  output=(char *)malloc(1024);

  sprintf(output,"%s.xz.nc",rootname);

  rmask=-9999.;
  status=create_ncfile_xy(output, grid.nx, grid.ny, grid,rmask);
  if(status !=0) goto error;

  status=nc_open(output,NC_WRITE,&ncid);
  if(status !=0) goto error;

/**----------------------------------------------------------------------------
  obsolete call, will be suppressed */
  status=nc_writeaxis(ncid,grid);
  if(status !=0) goto error;
   
  rbuffer=(float *)malloc(nndes*sizeof(float));
  rbuf=(float *)malloc(grid.nx*grid.ny*sizeof(float));

  count=0;
  do {
  count++;
  sprintf(sectionfile,"%s.%3.3d.s2r",rootname,count);
   printf("reading %s\n",sectionfile);
   //status=quoddy_loadr1(sectionfile, mesh.nvtxs, rbuffer,column); ERREUR avec archive.h Thierry
   status=quoddy_loadr1(sectionfile, mesh.nvtxs, rbuffer);
    if(status!=0) {
      printf("cannot read %s\n",sectionfile);
      goto finished;
      }
    status=fe_map(mesh,rbuffer,grid,elts,rbuf,rmask);
    if(status!=0) goto error;
/**----------------------------------------------------------------------------
    obsolete call, will be suppressed */
    status=nc_write_r1(ncid, grid, 3, count-1, rbuf,rmask);
    if(status!=0) goto error;
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
