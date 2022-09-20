

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
#include <errno.h>

#include "tools-structures.h"
#include "constants.h"
#include "netcdf-proto.h"

#include "sym-io.def"
#include "sym-io.h"
#include "geo.h"
#include "maths.h"
#include "map.h"


   /* variable ids */
static   int lon_id;
static   int lat_id;
static   int landmask_id;
static   int topo_id;
static   int cresolution_id, sresolution_id;


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int symio_save_data(const char *filename, const grid_t & grid, const char xloc)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   poc_var_t gridVar;
//     var=gridVar;
//     var.init("utransport",NC_DOUBLE,"sea_water_x_transport_"+SNS,"m2 s-1");
//     poc_put_vara(oP,var,fu,1);
//     var.init("vtransport",NC_DOUBLE,"sea_water_y_transport_"+SNS,"m2 s-1");
//     poc_put_vara(oP,var,fv,1);
//
//   status=poc_save_grid(filename, &gridVar, grid,1,0, xloc, 0.f, xloc, 0.f);
//   if(status)__NC_CHKERR_LINE_FILE__(status,"poc_save_grid() error");
//   status=poc_def_att(filename,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
//
//   return(status);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int symio_createfile(char *filename,size_t x_len, size_t y_len, const grid_t & grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/
{
   int  ncid;	/* netCDF id */

   /* dimension ids */
   int x_dim;
   int y_dim;


   /* variable shapes */
   int lon_dims[RANK_lon];
   int lat_dims[RANK_lat];
   int landmask_dims[RANK_landmask];
   int topo_dims[RANK_topo];
   int resolution_dims[RANK_resolution];

   /* attribute vectors */
   double lon_valid_min[1];
   double lon_valid_max[1];
   double lat_valid_min[1];
   double lat_valid_max[1];
   float landmask_missing_value[1];
   float landmask__FillValue[1];
   float landmask_scale_factor[1];
   float landmask_add_offset[1];
   float topo_missing_value[1];
   float topo__FillValue[1];
   float topo_scale_factor[1];
   float topo_add_offset[1];
   float resolution_missing_value[1];
   float resolution__FillValue[1];
   float resolution_scale_factor[1];
   float resolution_add_offset[1];

   /* enter define mode */
   int stat = nc_create(filename, NC_CLOBBER, &ncid);
   nc_check_error(stat,__LINE__,__FILE__,1);

   /* define dimensions */
   stat = nc_def_dim(ncid, "x", x_len, &x_dim);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_def_dim(ncid, "y", y_len, &y_dim);
   nc_check_error(stat,__LINE__,__FILE__,1);

   /* define variables */

   lon_dims[0] = y_dim;
   lon_dims[1] = x_dim;
   stat = nc_def_var(ncid, "lon", NC_DOUBLE, RANK_lon, lon_dims, &lon_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   lat_dims[0] = y_dim;
   lat_dims[1] = x_dim;
   stat = nc_def_var(ncid, "lat", NC_DOUBLE, RANK_lat, lat_dims, &lat_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   landmask_dims[0] = y_dim;
   landmask_dims[1] = x_dim;
   stat = nc_def_var(ncid, "landmask", NC_FLOAT, RANK_landmask, landmask_dims, &landmask_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   topo_dims[0] = y_dim;
   topo_dims[1] = x_dim;
   stat = nc_def_var(ncid, "topo", NC_FLOAT, RANK_topo, topo_dims, &topo_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   resolution_dims[0] = y_dim;
   resolution_dims[1] = x_dim;
   stat = nc_def_var(ncid, "c_resolution", NC_FLOAT, RANK_resolution, resolution_dims, &cresolution_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   resolution_dims[0] = y_dim;
   resolution_dims[1] = x_dim;
   stat = nc_def_var(ncid, "s_resolution", NC_FLOAT, RANK_resolution, resolution_dims, &sresolution_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   /* assign attributes */
   stat = nc_put_att_text(ncid, lon_id, "units", 12, "degrees_east");
   nc_check_error(stat,__LINE__,__FILE__,1);
   lon_valid_min[0] = -180;
   stat = nc_put_att_double(ncid, lon_id, "valid_min", NC_DOUBLE, 1, lon_valid_min);
   nc_check_error(stat,__LINE__,__FILE__,1);
   lon_valid_max[0] = 180;
   stat = nc_put_att_double(ncid, lon_id, "valid_max", NC_DOUBLE, 1, lon_valid_max);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, lon_id, "long_name", 9, "longitude");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, lon_id, "standard_name", 9, "longitude");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, lon_id, "nav_model", 12, "Default grid");
   nc_check_error(stat,__LINE__,__FILE__,1);

   stat = nc_put_att_text(ncid, lat_id, "units", 13, "degrees_north");
   nc_check_error(stat,__LINE__,__FILE__,1);
   lat_valid_min[0] = -90;
   stat = nc_put_att_double(ncid, lat_id, "valid_min", NC_DOUBLE, 1, lat_valid_min);
   nc_check_error(stat,__LINE__,__FILE__,1);
   lat_valid_max[0] = 90;
   stat = nc_put_att_double(ncid, lat_id, "valid_max", NC_DOUBLE, 1, lat_valid_max);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, lat_id, "long_name", 8, "latitude");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, lat_id, "standard_name", 8, "latitude");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, lat_id, "nav_model", 12, "Default grid");
   nc_check_error(stat,__LINE__,__FILE__,1);

   stat = nc_put_att_text(ncid, landmask_id, "units", 4, "none");
   nc_check_error(stat,__LINE__,__FILE__,1);
   landmask_missing_value[0] = 1e+35;
   stat = nc_put_att_float(ncid, landmask_id, "missing_value", NC_FLOAT, 1, landmask_missing_value);
   nc_check_error(stat,__LINE__,__FILE__,1);
   landmask__FillValue[0] = 1e+35;
   stat = nc_put_att_float(ncid, landmask_id, "_FillValue", NC_FLOAT, 1, landmask__FillValue);
   nc_check_error(stat,__LINE__,__FILE__,1);
   landmask_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, landmask_id, "scale_factor", NC_FLOAT, 1, landmask_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__,1);
   landmask_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, landmask_id, "add_offset", NC_FLOAT, 1, landmask_add_offset);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, landmask_id, "long_name", 9, "land mask");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, landmask_id, "standard_name", 9, "land_mask");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, landmask_id, "short_name", 8, "landmask");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, landmask_id, "content", 2, "YX");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, landmask_id, "associate", 7, "lat lon");
   nc_check_error(stat,__LINE__,__FILE__,1);

   stat = nc_put_att_text(ncid, topo_id, "units", 5, "meter");
   nc_check_error(stat,__LINE__,__FILE__,1);
   topo_missing_value[0] = 1e+35;
   stat = nc_put_att_float(ncid, topo_id, "missing_value", NC_FLOAT, 1, topo_missing_value);
   nc_check_error(stat,__LINE__,__FILE__,1);
   topo__FillValue[0] = 1e+35;
   stat = nc_put_att_float(ncid, topo_id, "_FillValue", NC_FLOAT, 1, topo__FillValue);
   nc_check_error(stat,__LINE__,__FILE__,1);
   topo_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, topo_id, "scale_factor", NC_FLOAT, 1, topo_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__,1);
   topo_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, topo_id, "add_offset", NC_FLOAT, 1, topo_add_offset);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, topo_id, "long_name", 10, "topography");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, topo_id, "standard_name", 10, "topography");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, topo_id, "short_name", 4, "topo");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, topo_id, "content", 2, "YX");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, topo_id, "associate", 7, "lat lon");
   nc_check_error(stat,__LINE__,__FILE__,1);

   stat = nc_put_att_text(ncid, cresolution_id, "units", 5, "meter");
   nc_check_error(stat,__LINE__,__FILE__,1);
   resolution_missing_value[0] = 1e+35;
   stat = nc_put_att_float(ncid, cresolution_id, "missing_value", NC_FLOAT, 1, resolution_missing_value);
   nc_check_error(stat,__LINE__,__FILE__,1);
   resolution__FillValue[0] = 1e+35;
   stat = nc_put_att_float(ncid, cresolution_id, "_FillValue", NC_FLOAT, 1, resolution__FillValue);
   nc_check_error(stat,__LINE__,__FILE__,1);
   resolution_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, cresolution_id, "scale_factor", NC_FLOAT, 1, resolution_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__,1);
   resolution_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, cresolution_id, "add_offset", NC_FLOAT, 1, resolution_add_offset);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cresolution_id, "long_name", 10, "resolution");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cresolution_id, "standard_name", 10, "resolution");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cresolution_id, "short_name", 4, "resolution");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cresolution_id, "content", 2, "YX");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cresolution_id, "associate", 7, "lat lon");
   nc_check_error(stat,__LINE__,__FILE__,1);

   stat = nc_put_att_text(ncid, sresolution_id, "units", 5, "meter");
   nc_check_error(stat,__LINE__,__FILE__,1);
   resolution_missing_value[0] = 1e+35;
   stat = nc_put_att_float(ncid, sresolution_id, "missing_value", NC_FLOAT, 1, resolution_missing_value);
   nc_check_error(stat,__LINE__,__FILE__,1);
   resolution__FillValue[0] = 1e+35;
   stat = nc_put_att_float(ncid, sresolution_id, "_FillValue", NC_FLOAT, 1, resolution__FillValue);
   nc_check_error(stat,__LINE__,__FILE__,1);
   resolution_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, sresolution_id, "scale_factor", NC_FLOAT, 1, resolution_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__,1);
   resolution_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, sresolution_id, "add_offset", NC_FLOAT, 1, resolution_add_offset);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, sresolution_id, "long_name", 10, "resolution");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, sresolution_id, "standard_name", 10, "resolution");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, sresolution_id, "short_name", 4, "resolution");
   nc_check_error(stat,__LINE__,__FILE__,1);
/* *------------------------------------------------------------------------
    attribute name changed to comply with CF standard*/
//   stat = nc_put_att_text(ncid, sresolution_id, "axis", 2, "YX");
   stat = nc_put_att_text(ncid, sresolution_id, "content", 2, "YX");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, sresolution_id, "associate", 7, "lat lon");
   nc_check_error(stat,__LINE__,__FILE__,1);

  stat = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", 21, "CF 1.0 revised by POC");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", 7, "regular");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "file_name", 12, "testhist1.nc");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, NC_GLOBAL, "production", 24, "merge produced this file");
   nc_check_error(stat,__LINE__,__FILE__,1);

   /* leave define mode */
   stat = nc_enddef (ncid);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_close(ncid);
   nc_check_error(stat,__LINE__,__FILE__,1);

   stat=nc_open(filename,NC_WRITE,&ncid);
   stat=nc_put_var_double(ncid,lon_id,grid.x);
   stat=nc_put_var_double(ncid,lat_id,grid.y);
   stat = nc_close(ncid);

   return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int symio_writefile(char *filename, float *landmask, float *topo, float *cresolution, float *sresolution)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
   int  ncid;

   int stat=nc_open(filename,NC_WRITE,&ncid);
   stat=nc_put_var_float(ncid,landmask_id,landmask);
   stat=nc_put_var_float(ncid,topo_id,topo);
   
   if(cresolution!=0) stat=nc_put_var_float(ncid,cresolution_id,cresolution);
   if(sresolution!=0) stat=nc_put_var_float(ncid,sresolution_id,sresolution);

   stat = nc_close(ncid);

   return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int savenodes(char *filename, serie_t metadata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in=NULL;
  int i,nndes;

  if ((in = fopen(filename, "w")) == NULL)
      return(-1);

  fprintf(in, "%d \n", metadata.count);
  fprintf(in, "%d \n", 0);
  fprintf(in, "%d \n", metadata.count);

  nndes=metadata.count;


  for (i=0; i<nndes; i++) {
    fprintf(in, "%f %f %f\n",metadata.data[i].lon,metadata.data[i].lat,metadata.data[i].values[0]);
    }
  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_rotate(geo_t *projection, const grid_t & cgrid, grid_t & sgrid, double pole_t, double pole_p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, status;
  matrix3x3_t R1,R2;
  vector3_t u,w;
  geo_t local;

  R1=math_rotation3D_init( 'y', (90.0-pole_p)*d2r);
  R2=math_rotation3D_init( 'z', pole_t*d2r);

  sgrid.dx=cgrid.dx;
  sgrid.dy=cgrid.dy;
  sgrid.nx=cgrid.nx;
  sgrid.ny=cgrid.ny;

  sgrid.modeH=2;

  sgrid.x=new double [sgrid.nx*sgrid.ny];
  sgrid.y=new double [sgrid.nx*sgrid.ny];
  
  local=geo_mercator_init(pole_t,pole_p,projection->radius);

  for(j=0;j<sgrid.ny;j++) {
    for(i=0;i<sgrid.nx;i++) {
      k=sgrid.nx*j+i;
      status=geo_mercator_inverse(*projection,&(sgrid.x[k]),&(sgrid.y[k]),cgrid.x[k],cgrid.y[k]);
      status=math_rotation3Ddeg('D', R1, R2, &(sgrid.x[k]),&(sgrid.y[k]));
      }
    }

  *projection=local;
  
  range_t<double> range;
  
  range=range_t<double>(sgrid.x,sgrid.Hsize());
  sgrid.xmin=range.min;
  sgrid.xmax=range.max;
  range=range_t<double>(sgrid.y,sgrid.Hsize());
  sgrid.ymin=range.min;
  sgrid.ymax=range.max;

  return(0);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

grid_t get_z2ugrid(grid_t zgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,k1,k2;
  double dx;
  grid_t ugrid;

  ugrid.dx=zgrid.dx;
  ugrid.dy=zgrid.dy;
  ugrid.nx=zgrid.nx;
  ugrid.ny=zgrid.ny;
  ugrid.nz=zgrid.nz;
  ugrid.modeH=0;

  ugrid.x=new double[ugrid.nx*ugrid.ny];
  ugrid.y=new double[ugrid.nx*ugrid.ny];

  ugrid.ymin= +1.e+35;
  ugrid.ymax= -1.e+35;
  ugrid.xmin= +1.e+35;
  ugrid.xmax= -1.e+35;

  for(j=0;j<ugrid.ny;j++) {
    for(i=0;i<ugrid.nx;i++) {
      k=ugrid.nx*j+i;
      if(i==0) {
        k1=zgrid.nx*j+i;
        k2=zgrid.nx*j+i+1;
        }
      else {
        k1=zgrid.nx*j+i-1;
        k2=zgrid.nx*j+i;
        }
      dx=zgrid.x[k2]-zgrid.x[k1];
      ugrid.x[k]=zgrid.x[k]-0.5*dx;
      ugrid.y[k]=zgrid.y[k];
      ugrid.ymin=min(ugrid.ymin,ugrid.y[k]);
      ugrid.ymax=max(ugrid.ymax,ugrid.y[k]);
      ugrid.xmin=min(ugrid.xmin,ugrid.x[k]);
      ugrid.xmax=max(ugrid.xmax,ugrid.x[k]);
      }
    }

  return(ugrid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

grid_t get_z2vgrid(grid_t zgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,k1,k2;
  double dy;
  grid_t vgrid;

  vgrid.dx=zgrid.dx;
  vgrid.dy=zgrid.dy;
  vgrid.nx=zgrid.nx;
  vgrid.ny=zgrid.ny;
  vgrid.nz=zgrid.nz;
  vgrid.modeH=0;

  vgrid.x=new double[vgrid.nx*vgrid.ny];
  vgrid.y=new double[vgrid.nx*vgrid.ny];

  vgrid.ymin= +1.e+35;
  vgrid.ymax= -1.e+35;
  vgrid.xmin= +1.e+35;
  vgrid.xmax= -1.e+35;

  for(j=0;j<vgrid.ny;j++) {
    for(i=0;i<vgrid.nx;i++) {
      k=vgrid.nx*j+i;
      if(j==0) {
        k1=zgrid.nx*j+i;
        k2=zgrid.nx*(j+1)+i;
        }
      else {
        k1=zgrid.nx*(j-1)+i;
        k2=zgrid.nx*j+i;
        }
      vgrid.x[k]=zgrid.x[k];
      dy=zgrid.y[k2]-zgrid.y[k1];
      vgrid.y[k]=zgrid.y[k]-0.5*dy;
      vgrid.ymin=min(vgrid.ymin,vgrid.y[k]);
      vgrid.ymax=max(vgrid.ymax,vgrid.y[k]);
      vgrid.xmin=min(vgrid.xmin,vgrid.x[k]);
      vgrid.xmax=max(vgrid.xmax,vgrid.x[k]);
      }
    }

  return(vgrid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

grid_t get_f2zgrid(grid_t fgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k;
  grid_t zgrid;

//   zgrid.dx=fgrid.dx;
//   zgrid.dy=fgrid.dy;
  zgrid.nx=fgrid.nx-1;
  zgrid.ny=fgrid.ny-1;
  zgrid.nz=fgrid.nz;
  zgrid.modeH=fgrid.modeH;

  zgrid.x=new double[zgrid.Hsize()];
  zgrid.y=new double[zgrid.Hsize()];

  zgrid.ymin= +1.e+35;
  zgrid.ymax= -1.e+35;
  zgrid.xmin= +1.e+35;
  zgrid.xmax= -1.e+35;

  for(j=0;j<zgrid.ny;j++) {
    for(i=0;i<zgrid.nx;i++) {
      k=zgrid.nx*j+i;
      int m1=fgrid.Hindex(i,j);
      int m2=fgrid.Hindex(i+1,j);
      int m3=fgrid.Hindex(i+1,j+1);
      int m4=fgrid.Hindex(i,j+1);
      zgrid.x[k]=0.25*(fgrid.x[m1]+fgrid.x[m2]+fgrid.x[m3]+fgrid.x[m4]);
      zgrid.y[k]=0.25*(fgrid.x[m1]+fgrid.x[m2]+fgrid.x[m3]+fgrid.x[m4]);
      zgrid.ymin=min(zgrid.ymin,zgrid.y[k]);
      zgrid.ymax=max(zgrid.ymax,zgrid.y[k]);
      zgrid.xmin=min(zgrid.xmin,zgrid.x[k]);
      zgrid.xmax=max(zgrid.xmax,zgrid.x[k]);
      }
    }

  return(zgrid);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int read_notebook(const char *input, ntbk_grid_t *notebook)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    n=1024;
  int    nitems;
  FILE   *in=NULL;
  char   line[1024];
  char   *test=NULL;

  in=fopen(input,"r");

  if(in==NULL) {
    printf("unable to open notebook file: %s, abort...\n",input);
    return(-1);
    }

  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%lf",&notebook->dXb);
  if(nitems<1){fclose(in);__TRAP_ERR_RETURN__(-1,1,"error reading %s\n",input);}
  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%lf",&notebook->dYb);
/*------------------------------------------------------------------------
  reference latitude for geographical projection */
  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%lf",&notebook->Phi_0);

/*------------------------------------------------------------------------
  longitude,latitude of I0,J0 */
  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%lf",&notebook->Longi_0);
  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%lf",&notebook->Latit_0);

  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%lf",&notebook->Angle_0);

  fgets(line,n,in);
  nitems=sscanf(line,"%lf",&notebook->I0);
  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%lf",&notebook->J0);

  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%lf",&notebook->RayonTerre);

  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%d",&notebook->TypeGrid);
  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%d",&notebook->I1D);

  line[0]='\0';
  fgets(line,n,in);
  nitems=sscanf(line,"%d %d %d",&notebook->MECO,&notebook->NECO,&notebook->NR);

/* *------------------------------------------------------------------------
  2010 add-on for grid distorsion*/
  test=fgets(line,n,in);
  if(test==0) {
    notebook->Pole_LON=0.0;
    }
  else {
    nitems=sscanf(line,"%lf",&notebook->Pole_LON);
    if(nitems==0) notebook->Pole_LON=0.0;
    }

  test=fgets(line,n,in);
  if(test==0) {
    notebook->Pole_LAT=90.0;
    }
  else {
    nitems=sscanf(line,"%lf",&notebook->Pole_LAT);
    if(nitems==0) notebook->Pole_LAT=90.0;
    }

  fclose(in);

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool is_integer(double value)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  bool chk;
  
  double tmp=floor(value+0.5);
  
  if(tmp==value) chk=true;
  else chk=false;
  
  return(chk);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int save_notebook(const char *input, const ntbk_grid_t & notebook)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    nitems;
  FILE   *in;

  in=fopen(input,"w");

  if(in==NULL) {
    int status=errno;
    printf("unable to open notebook file: %s, abort...\n",input);
    return status;
    }
  nitems=fprintf(in,"%lf %s\n",notebook.dXb,"DXB        Along Oi axis cellboxe size (m)");
  nitems=fprintf(in,"%lf %s\n",notebook.dYb,"DYB        Along Oj axis cellboxe size (m)");
/*------------------------------------------------------------------------
  reference latitude for geographical projection */
  nitems=fprintf(in,"%lf %s\n",notebook.Phi_0,"PHI0     Reference latitude for the Mercator projection");

/*------------------------------------------------------------------------
  longitude,latitude of I0,J0 */
  nitems=fprintf(in,"%lf %s\n",notebook.Longi_0,"LONGI0     longitude of grid point (I0,J0)");
  nitems=fprintf(in,"%lf %s\n",notebook.Latit_0,"LATIT0     latitude of grid point (I0,J0)");

  nitems=fprintf(in,"%lf %s\n",notebook.Angle_0,"not relevant in S25");

// /*------------------------------------------------------------------------------
//   testing */
//   fclose(in);
//   return(0);

  if(is_integer(notebook.I0))
    nitems=fprintf(in,"%d %s\n",nint(notebook.I0),"I0");
  else
    nitems=fprintf(in,"%lf %s\n",notebook.I0,"I0");

  if(is_integer(notebook.J0))
    nitems=fprintf(in,"%d %s\n",nint(notebook.J0),"J0");
  else
    nitems=fprintf(in,"%lf %s\n",notebook.J0,"J0");

  nitems=fprintf(in,"%lf %s\n",notebook.RayonTerre,"RAYONTERRE Earth radius (m)");

  nitems=fprintf(in,"%d %s\n",notebook.TypeGrid,"TYPEGRID   1=Mercator");
  nitems=fprintf(in,"%d %s\n",notebook.I1D,     "I1D        1=3D model");

  nitems=fprintf(in,"%d %d %d %s\n",notebook.MECO,notebook.NECO,notebook.NR,"iglb jglb kmax (as in module_parameter.F90)");

/**------------------------------------------------------------------------
  2010 add-on for grid distorsion*/
  nitems=fprintf(in,"%lf %s\n",notebook.Pole_LON,"POLE_LON   ! longitude (° decimal) of the grid pole");

  nitems=fprintf(in,"%lf %s\n",notebook.Pole_LAT,"POLE_LAT   ! latitude  (° decimal) of the grid pole");

  const char *comment="\n \
  ..... \n \
Note: a 360° grid is obtained with\n \
DXB=DYB=rayonterre*2.*pi/(imax-2)*cos(latit0(radians))\n \
and the special mpi option in module_parameter.F90\n";

  nitems=fprintf(in,"%s\n",comment);

  fclose(in);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int process_cartesian_notebook(const ntbk_grid_t & notebook, grid_t *cgrid, grid_t & sgrid, geo_t *projection, bool extend /*, bool force_consistency*/)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  build cartesian and geocentric grid from cartesian specification
  
  use non-standard mercator projection (might be an issue)
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    i,j,k;
  double I0,J0;
  int    status;
  double lon,lat,angle,radius,x,y;
  double tx,ty,nx,ny;

  cgrid->dx=notebook.dXb;
  cgrid->dy=notebook.dYb;

  projection->lat=notebook.Phi_0;

  lon=notebook.Longi_0;
  lat=notebook.Latit_0;

  angle=notebook.Angle_0;

  I0=notebook.I0;
  J0=notebook.J0;

  radius=notebook.RayonTerre;

  cgrid->nx=notebook.MECO;
  cgrid->ny=notebook.NECO;
  cgrid->nz=notebook.NR;

/* *------------------------------------------------------------------------
  [1:M,1:N] are extended to [0:M+1,0:N+1]  */
  if(extend) {
    cgrid->nx+=2;
    cgrid->ny+=2;
    I0++;
    J0++;
    }
/*------------------------------------------------------------------------
  initialize mercator projection */
  *projection=geo_mercator_init(lon,projection->lat,radius);

/*------------------------------------------------------------------------
  cartesian coordinates of I0,J0 */
  status=geo_mercator_directe(*projection,lon,lat,&x,&y);

/*------------------------------------------------------------------------
  cartesian grid orientation and axis unitary vectors */
  angle*=-d2r;

  tx=cos(angle);
  ty=sin(angle);

  nx=-ty;
  ny= tx;

  cgrid->x=new double[cgrid->nx*cgrid->ny];
  cgrid->y=new double[cgrid->nx*cgrid->ny];
  cgrid->z=NULL;

  cgrid->modeH=2;

  cgrid->xmin=x-I0*cgrid->dx*tx-J0*cgrid->dy*nx;
  cgrid->ymin=y-I0*cgrid->dx*ty-J0*cgrid->dy*ny;

  for(j=0;j<cgrid->ny;j++) {
    for(i=0;i<cgrid->nx;i++) {
      k=cgrid->nx*j+i;
      cgrid->x[k]= cgrid->xmin+i*tx*cgrid->dx+j*nx*cgrid->dy;
      cgrid->y[k]= cgrid->ymin+i*ty*cgrid->dx+j*ny*cgrid->dy;
      }
    }

  k=cgrid->nx*cgrid->ny-1;
  cgrid->xmax=cgrid->x[k];
  cgrid->ymax=cgrid->y[k];

  if((notebook.Pole_LON!=0) || (notebook.Pole_LAT!=90.)) {
    status= map_rotate(projection, *cgrid, sgrid, notebook.Pole_LON, notebook.Pole_LAT);
    }
  else {
    sgrid=map_get_spherical(*projection,*cgrid);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  implemented for mesh-academic
  
  commented
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   char parameters[1024];
//   projPJ pj=assign_projection(0, &projection->lat, &lon, parameters);
//   cgrid->projection=pj;
//   cgrid->proj4options=poc_strdup(parameters);
//   
// //   sgrid=map_get_spherical(pj,*cgrid);
//   for(int m=0; m< cgrid->Hsize(); m++) {
//     cgrid->x[m]=sgrid.x[m];
//     cgrid->y[m]=sgrid.y[m];
//     }
//   status=geo_to_projection (cgrid->proj4options, cgrid->x, cgrid->y, cgrid->Hsize());
//   
//   status=map_minmax(cgrid);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int process_periodic_notebook(const ntbk_grid_t & notebook, grid_t *cgrid, grid_t & sgrid, geo_t *projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k;
  double I0,J0;
  int    status;
  double lon,lat,angle,radius,x,y;
  double tx,ty,nx,ny;
  double dxb,dyb;

  projection->lat=notebook.Phi_0;

  lon=notebook.Longi_0;
  lat=notebook.Latit_0;

  angle=notebook.Angle_0;

  I0=notebook.I0;
  J0=notebook.J0;

  radius=notebook.RayonTerre;

  cgrid->nx=notebook.MECO;
  cgrid->ny=notebook.NECO;
  cgrid->nz=notebook.NR;

  dxb=cos(lat*d2r)*radius*2.*M_PI/double(cgrid->nx-2);
  dyb=dxb;

  cgrid->dx=dxb;
  cgrid->dy=dyb;

/* *------------------------------------------------------------------------
  [1:M,1:N] are extended to [0:M+1,0:N+1]  */
  cgrid->nx+=2;
  cgrid->ny+=2;

/*------------------------------------------------------------------------
  initialize mercator projection */
  *projection=geo_mercator_init(lon,projection->lat,radius);

/*------------------------------------------------------------------------
  cartesian coordinates of I0,J0 */
  status=geo_mercator_directe(*projection,lon,lat,&x,&y);

/*------------------------------------------------------------------------
  cartesian grid orientation and axis unitary vectors */
  angle*=-d2r;

  tx=cos(angle);
  ty=sin(angle);

  nx=-ty;
  ny= tx;

  cgrid->x=new double[cgrid->nx*cgrid->ny];
  cgrid->y=new double[cgrid->nx*cgrid->ny];
  cgrid->z=NULL;

  cgrid->modeH=2;

  cgrid->xmin=x-I0*cgrid->dx*tx-J0*cgrid->dy*nx;
  cgrid->ymin=y-I0*cgrid->dx*ty-J0*cgrid->dy*ny;

  for(j=0;j<cgrid->ny;j++) {
    for(i=0;i<cgrid->nx;i++) {
      k=cgrid->nx*j+i;
      cgrid->x[k]= cgrid->xmin+i*tx*cgrid->dx+j*nx*cgrid->dy;
      cgrid->y[k]= cgrid->ymin+i*ty*cgrid->dx+j*ny*cgrid->dy;
      }
    }

  k=cgrid->nx*cgrid->ny-1;
  cgrid->xmax=cgrid->x[k];
  cgrid->ymax=cgrid->y[k];

  if((notebook.Pole_LON!=0) || (notebook.Pole_LAT!=90.)) {
    status= map_rotate(projection, *cgrid, sgrid, notebook.Pole_LON, notebook.Pole_LAT);
    }
  else {
    sgrid=map_get_spherical(*projection,*cgrid);
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int process_spherical_notebook(const ntbk_grid_t & notebook, grid_t *cgrid, geo_t *projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k;
  double I0,J0;
  double lon,lat,angle,radius,x,y;
  double tx,ty,nx,ny;


  cgrid->dx=notebook.dXb;
  cgrid->dy=notebook.dYb;

  lon=notebook.Longi_0;
  lat=notebook.Latit_0;

  angle=notebook.Angle_0;

  I0=notebook.I0;
  J0=notebook.J0;

  radius=notebook.RayonTerre;

  cgrid->nx=notebook.MECO;
  cgrid->ny=notebook.NECO;
  cgrid->nz=notebook.NR;

/* *------------------------------------------------------------------------
  [1:M,1:N] are extended to [0:M+1,0:N+1]  */
  cgrid->nx+=2;
  cgrid->ny+=2;

/*------------------------------------------------------------------------
  initialize mercator projection */
  *projection=geo_mercator_init(lon,projection->lat,radius);

  x=lon;
  y=lat;

/*------------------------------------------------------------------------
  cartesian grid orientation and axis unitary vectors */
  angle*=-d2r;

  tx=cos(angle);
  ty=sin(angle);

  nx=-ty;
  ny= tx;

  cgrid->x=new double[cgrid->nx*cgrid->ny];
  cgrid->y=new double[cgrid->nx*cgrid->ny];
  cgrid->z=NULL;

  cgrid->modeH=2;

  cgrid->xmin=x-I0*cgrid->dx*tx-J0*cgrid->dy*nx;
  cgrid->ymin=y-I0*cgrid->dx*ty-J0*cgrid->dy*ny;

  for(j=0;j<cgrid->ny;j++) {
    for(i=0;i<cgrid->nx;i++) {
      k=cgrid->nx*j+i;
      cgrid->x[k]= cgrid->xmin+i*tx*cgrid->dx+j*nx*cgrid->dy;
      cgrid->y[k]= cgrid->ymin+i*ty*cgrid->dx+j*ny*cgrid->dy;
      }
    }

  k=cgrid->nx*cgrid->ny-1;
  cgrid->xmax=cgrid->x[k];
  cgrid->ymax=cgrid->y[k];

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int load_notebook(const char *input, grid_t *cgrid, grid_t *sgrid, geo_t *projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  get grid from notebook file

  see http://sirocco.omp.obs-mip.fr/outils/Symphonie/Documentation/SymphonieDocNotebook.htm#grid

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int    status;
  ntbk_grid_t notebook;

  status=read_notebook(input, &notebook);
  if(status!=0) return(status);
  
  bool extend=false;
  status=notebook2grids(notebook, cgrid,sgrid, extend, projection);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int notebook2grids(const ntbk_grid_t & notebook, grid_t *cgrid, grid_t *sgrid, bool extend, geo_t *projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
   build grid from notebook instructions
   
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  
  geo_t proj_;
  if(projection==0) projection=&proj_;

  switch (notebook.TypeGrid) {
/*------------------------------------------------------------------------------
    Cartesian, possibly polar, grid definition */
    case 1:
      status=process_cartesian_notebook(notebook, cgrid, *sgrid, projection, extend);
//      *sgrid=map_get_spherical(*projection,*cgrid);
      break;

/*------------------------------------------------------------------------------
    Geocentric grid definition */
    case 2:
      status=process_spherical_notebook(notebook, sgrid, projection);
      *cgrid=map_get_cartesian(*projection,*sgrid);
      break;

/*------------------------------------------------------------------------------
    Periodic grid definition */
    case 11:
      status=process_periodic_notebook(notebook, cgrid, *sgrid, projection);
      *sgrid=map_get_spherical(*projection,*cgrid);
      break;

    default:
      status=-1;
      break;
    }

  cgrid->connex=1;
  sgrid->connex=1;

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int symphonie_loadmask(char *input,grid_t sgrid,float *landmask,float  *topo)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l;
  FILE *in=NULL;
  double x,y;

  in=fopen(input,"r");
  if(in==NULL) return(-1);

  for(i=1;i<sgrid.nx-1;i++) {
    for(j=1;j<sgrid.ny-1;j++) {
      k=sgrid.nx*j+i;
      if (j%400==0) fprintf(in,"\n");
      fscanf(in,"%d %d %f %lf %1f",&k,&l,&landmask[k],&x,&y);
      }
    }

  fclose(in);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int symphonie_savemask(char *output,grid_t sgrid,float *landmask,float  *topo)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k;
  int flag;
  FILE *out=NULL;

  out=fopen(output,"w");
  if(out==NULL) return(-1);

  for(i=1;i<sgrid.nx-1;i++) {
    for(j=1;j<sgrid.ny-1;j++) {
      k=sgrid.nx*j+i;
      if (j%400==0) fprintf(out,"\n");
      switch ((int) landmask[k]) {
        case 0:
          flag=0;
          break;
        case -1:
          flag=0;
          break;
        case 1:
          flag=1;
          break;
        }
      fprintf(out,"%1d",flag);
      }
    fprintf(out,"\n");
    }
  for(i=1;i<sgrid.nx-1;i++) {
    for(j=1;j<sgrid.ny-1;j++) {
      k=sgrid.nx*j+i;
      if (j%400==0) fprintf(out,"\n");
      fprintf(out," %6.1f",topo[k]);
      }
    fprintf(out,"\n");
    }

  fclose(out);
  return(0);
}
