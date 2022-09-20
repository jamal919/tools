#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"
 
#include "poc-netcdf.def"
#  define RANK_mss 2
#  define RANK_cls 2
#  define RANK_innovation 2
#  define RANK_geoid 2
#  define RANK_topo 2
  /* IB = mog2D medsea */
#include "netcdf-proto.h"

#define MOG 1
#define TIDE 18 /* regional model */
#include "poc-time.h"
#include "tides.h"
#include "tides.def"
#include "filter.h"
#include "statistic.h"
#include "polygones.h"
#include "geo.h"
#include "map.h"
#include "rutin.h"

#include "mss.v2.h"

   //  variable ids
   int lon_id;
   int lat_id;
   int mss_id;
   int cls_id;
   int innovation_id;
   int geoid_id;
   int topo_id;

/*------------------------------------------------------------------------
  
  Nom         :  createfile

  --------------------------------------------------------------------------*/

int mss_createfile(char *filename,size_t x_len, size_t y_len, grid_t grid)
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/

{		
   int  ncid;			/* netCDF id */

   /* dimension ids */
   int x_dim;
   int y_dim;



   /* variable shapes */
   int lon_dims[RANK_lon];
   int lat_dims[RANK_lat];
   int mss_dims[RANK_mss];
   int cls_dims[RANK_cls];
   int innovation_dims[RANK_innovation];
   int geoid_dims[RANK_geoid];
   int topo_dims[RANK_topo];

   /* attribute vectors */
   double lon_valid_min[1];
   double lon_valid_max[1];
   double lat_valid_min[1];
   double lat_valid_max[1];
   float mss_missing_value[1];
   float mss__FillValue[1];
   float mss_scale_factor[1];
   float mss_add_offset[1];
   float cls_missing_value[1];
   float cls__FillValue[1];
   float cls_scale_factor[1];
   float cls_add_offset[1];
   float innovation_missing_value[1];
   float innovation__FillValue[1];
   float innovation_scale_factor[1];
   float innovation_add_offset[1];
   float geoid_missing_value[1];
   float geoid__FillValue[1];
   float geoid_scale_factor[1];
   float geoid_add_offset[1];
   float topo_missing_value[1];
   float topo__FillValue[1];
   float topo_scale_factor[1];
   float topo_add_offset[1];

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

   mss_dims[0] = y_dim;
   mss_dims[1] = x_dim;
   stat = nc_def_var(ncid, "mss", NC_FLOAT, RANK_mss, mss_dims, &mss_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   cls_dims[0] = y_dim;
   cls_dims[1] = x_dim;
   stat = nc_def_var(ncid, "cls", NC_FLOAT, RANK_cls, cls_dims, &cls_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   innovation_dims[0] = y_dim;
   innovation_dims[1] = x_dim;
   stat = nc_def_var(ncid, "innovation", NC_FLOAT, RANK_innovation, innovation_dims, &innovation_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   geoid_dims[0] = y_dim;
   geoid_dims[1] = x_dim;
   stat = nc_def_var(ncid, "geoid", NC_FLOAT, RANK_geoid, geoid_dims, &geoid_id);
   nc_check_error(stat,__LINE__,__FILE__,1);

   topo_dims[0] = y_dim;
   topo_dims[1] = x_dim;
   stat = nc_def_var(ncid, "topo", NC_FLOAT, RANK_topo, topo_dims, &topo_id);
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
   stat = nc_put_att_text(ncid, mss_id, "units", 5, "meter");
   nc_check_error(stat,__LINE__,__FILE__,1);
   mss_missing_value[0] = 1e+35;
   stat = nc_put_att_float(ncid, mss_id, "missing_value", NC_FLOAT, 1, mss_missing_value);
   nc_check_error(stat,__LINE__,__FILE__,1);
   mss__FillValue[0] = 1e+35;
   stat = nc_put_att_float(ncid, mss_id, "_FillValue", NC_FLOAT, 1, mss__FillValue);
   nc_check_error(stat,__LINE__,__FILE__,1);
   mss_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, mss_id, "scale_factor", NC_FLOAT, 1, mss_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__,1);
   mss_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, mss_id, "add_offset", NC_FLOAT, 1, mss_add_offset);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, mss_id, "long_name", 14, "mean sea level");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, mss_id, "standard_name", 14, "mean_sea_level");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, mss_id, "short_name", 3, "mss");
   nc_check_error(stat,__LINE__,__FILE__,1);

   stat = nc_put_att_text(ncid, cls_id, "units", 5, "meter");
   nc_check_error(stat,__LINE__,__FILE__,1);
   cls_missing_value[0] = 1e+35;
   stat = nc_put_att_float(ncid, cls_id, "missing_value", NC_FLOAT, 1, cls_missing_value);
   nc_check_error(stat,__LINE__,__FILE__,1);
   cls__FillValue[0] = 1e+35;
   stat = nc_put_att_float(ncid, cls_id, "_FillValue", NC_FLOAT, 1, cls__FillValue);
   nc_check_error(stat,__LINE__,__FILE__,1);
   cls_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, cls_id, "scale_factor", NC_FLOAT, 1, cls_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__,1);
   cls_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, cls_id, "add_offset", NC_FLOAT, 1, cls_add_offset);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cls_id, "long_name", 18, "cls mean sea level");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cls_id, "standard_name", 18, "cls_mean_sea_level");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cls_id, "short_name", 8, "cls_msss");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cls_id, "content", 2, "YX");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, cls_id, "associate", 7, "lat lon");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, innovation_id, "units", 5, "meter");
   nc_check_error(stat,__LINE__,__FILE__,1);
   innovation_missing_value[0] = 1e+35;
   stat = nc_put_att_float(ncid, innovation_id, "missing_value", NC_FLOAT, 1, innovation_missing_value);
   nc_check_error(stat,__LINE__,__FILE__,1);
   innovation__FillValue[0] = 1e+35;
   stat = nc_put_att_float(ncid, innovation_id, "_FillValue", NC_FLOAT, 1, innovation__FillValue);
   nc_check_error(stat,__LINE__,__FILE__,1);
   innovation_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, innovation_id, "scale_factor", NC_FLOAT, 1, innovation_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__,1);
   innovation_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, innovation_id, "add_offset", NC_FLOAT, 1, innovation_add_offset);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, innovation_id, "long_name", 25, "mean sea level innovation");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, innovation_id, "standard_name", 25, "mean_sea_level_innovation");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, innovation_id, "short_name", 10, "innovation");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, innovation_id, "content", 2, "YX");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, innovation_id, "associate", 7, "lat lon");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, geoid_id, "units", 5, "meter");
   nc_check_error(stat,__LINE__,__FILE__,1);
   geoid_missing_value[0] = 1e+35;
   stat = nc_put_att_float(ncid, geoid_id, "missing_value", NC_FLOAT, 1, geoid_missing_value);
   nc_check_error(stat,__LINE__,__FILE__,1);
   geoid__FillValue[0] = 1e+35;
   stat = nc_put_att_float(ncid, geoid_id, "_FillValue", NC_FLOAT, 1, geoid__FillValue);
   nc_check_error(stat,__LINE__,__FILE__,1);
   geoid_scale_factor[0] = 1;
   stat = nc_put_att_float(ncid, geoid_id, "scale_factor", NC_FLOAT, 1, geoid_scale_factor);
   nc_check_error(stat,__LINE__,__FILE__,1);
   geoid_add_offset[0] = 0;
   stat = nc_put_att_float(ncid, geoid_id, "add_offset", NC_FLOAT, 1, geoid_add_offset);
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, geoid_id, "long_name", 5, "geoid");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, geoid_id, "standard_name", 5, "geoid");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, geoid_id, "short_name", 5, "geoid");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, geoid_id, "content", 2, "YX");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, geoid_id, "associate", 7, "lat lon");
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
   stat = nc_put_att_text(ncid, topo_id, "long_name", 17, "bottom topography");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, topo_id, "standard_name", 17, "bottom_topography");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, topo_id, "short_name", 4, "topo");
   nc_check_error(stat,__LINE__,__FILE__,1);


   stat = nc_put_att_text(ncid, mss_id, "content", 2, "YX");
   nc_check_error(stat,__LINE__,__FILE__,1);
   stat = nc_put_att_text(ncid, mss_id, "associate", 7, "lat lon");
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


/*------------------------------------------------------------------------
  
  Nom         :  writefile

  --------------------------------------------------------------------------*/

int writefile(char *filename, float *mss, float *cls, float *ino)

{		
   int  ncid;			/* netCDF id */

 
   int stat=nc_open(filename,NC_WRITE,&ncid);
   stat=nc_put_var_float(ncid,mss_id,mss);
   stat=nc_put_var_float(ncid,cls_id,cls);
   stat=nc_put_var_float(ncid,innovation_id,ino);
   stat = nc_close(ncid);

   return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int savenodes(char *filename, serie_t metadata)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  FILE *in;
  float dummy;
  vertex_t *set;
  int i,j,nndes,maxn,ndum,ninteriors,nexteriors;

  if ((in = fopen(filename, "w")) == NULL)
      return(-1);

  fprintf(in, "%d \n", metadata.count);
  fprintf(in, "%d \n", 0);
  fprintf(in, "%d \n", metadata.count);

  nndes=metadata.count;


  for (i=0; i<nndes; i++) {
      if(metadata.data[i].lon>180)metadata.data[i].lon-=360.0;
      fprintf(in, "%f %f %f\n",metadata.data[i].lon,metadata.data[i].lat,metadata.data[i].values[0]);
    }
  fclose(in);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

serie_t load_metadata_ref(char *input, int point)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,l,n;
  int    npoints,nitems,nmes,count,status;
  int    *keep,threshold;

  float  dum;
  float  mss,hf;

  double time,lon,lat;
  double *serie,*z,*t,*residuals,mean;

  char   line[1024],nodefile[1024],dummy[1024];

  FILE   *in,*out;

  data_t  *buffer;
  serie_t metadata;
  spectrum_t s;
  statistic_t sdum;

  metadata.count=0;
  metadata.data=NULL;

  in=fopen(input,"r");
  if(in==NULL) return(metadata);

  n=999;
  for(i=0;i<29;i++) fgets(line,n,in);
  nitems=fscanf(in,"%d\n",&npoints);

 read:
  for(i=0;i<1;i++) fgets(line,n,in);

  for(i=0;i<1;i++) fgets(line,n,in);
  nitems=sscanf(strstr(line," : "),"%s %d",dummy,&count);
  for(i=0;i<1;i++) fgets(line,n,in);
  nitems=sscanf(strstr(line," : "),"%s %lf",dummy,&lon);
  for(i=0;i<1;i++) fgets(line,n,in);
  nitems=sscanf(strstr(line," : "),"%s %lf",dummy,&lat);
  for(i=0;i<1;i++) fgets(line,n,in);
  nitems=sscanf(strstr(line," : "),"%s %d",dummy,&nmes);
 
  printf("file %s : %d %d %lf %lf\n",input,count,nmes,lon,lat);

  if(count!=point) {
      for(i=0;i<nmes;i++) fgets(line,n,in);
      goto read;
    }


  buffer=new data_t[nmes];
  k=0;
  while (k<nmes)
    {
      nitems=fscanf(in,"%lf",&(buffer[k].time));
      buffer[k].lon=lon;
      buffer[k].lat=lat;
      if(nitems!=1) break;
      
      for(j=0;j<24;j++) {
          nitems=fscanf(in,"%f",&(buffer[k].values[j]));
          if(nitems!=1) break;
        }
      k++;
    }
  fclose(in);

  metadata.count=nmes;
  metadata.mask=99.9999;
  metadata.data=buffer;

  keep=(int *)calloc(nmes,sizeof(int));
  serie=(double *)malloc(nmes*sizeof(double));
  t=(double *)malloc(nmes*sizeof(double));

  for(k=0;k<nmes;k++) {
      buffer[k].values[29]=metadata.mask;
      if(buffer[k].values[0] != metadata.mask){
        hf=buffer[k].values[MOG]+buffer[k].values[21];
        buffer[k].values[29]=buffer[k].values[0]-hf;
        buffer[k].values[28]=buffer[k].values[0]-buffer[k].values[MOG];
        }
    }

  printf("%s starts %s, finishes %s \n",input,sgetcnesdate(buffer[0].time*24.),sgetcnesdate(buffer[nmes-1].time*24.));
  sprintf(nodefile,"%s.nod",input);


/*-----------------------------------------------------------------------------
  remove annual signal*/

  status=spectrum_init(&s);
  status=addwave(&s, wSa);
  spectrum_terminate(s);

  residuals=(double *)malloc(nmes*sizeof(double));

  for(k=0,l=0;k<nmes;k++)
    if(buffer[k].values[0] != metadata.mask) {
      serie[l]=buffer[k].values[29];
      t[l]=buffer[k].time;
      l++;
      }
  count=l;
  status=harmonic_analysis_old(serie,t,residuals,&mean,count, s);

  for(k=0,l=0;k<nmes;k++)
    if(buffer[k].values[0] != metadata.mask) {
      buffer[k].values[29]=residuals[l];
      l++;
      }
  sdum=get_statistics(residuals, (double) metadata.mask, nmes);

  out=fopen("gauge.dat","w");
  for(k=0;k<nmes;k++)
    if(buffer[k].values[0] != metadata.mask)
      fprintf(out,"%lf %9.3f %9.3f\n",buffer[k].time,buffer[k].values[0]-mean,buffer[k].values[29]);
  fclose(out);

  free(serie); serie=NULL;
  free(t);     t=NULL;
  free(residuals); residuals=NULL;

  return(metadata);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int decale_metadata(data_t *data,int *keep,int nmes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n;


  for(j=0,n=0;j<nmes;j++) {
      if (keep[j]==1) continue;
      else
        {
          data[n].lon=data[j].lon;
          data[n].lat=data[j].lat;
          data[n].cycle=data[j].cycle;
          data[n].time=data[j].time;
          for(i=0;i<30;i++) data[n].values[i]=data[j].values[i];
          n++;
        }
    }

  return(n);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

serie_t load_metadata_raw(char *input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    i,j,k,n,npoints=0,nitems,nmes,cycle,status;
  int    interior,threshold,*keep;
  float  dum,tides,mss;
  double dt,cut;
  double time,lon,lat,*serie,*t,*residuals,mean;
  double *z;
  char   line[1024],nodefile[1024];
  FILE   *in;
  data_t  *buffer;
  serie_t metadata;
  spectrum_t s;
  statistic_t sdum;

  if( (in=fopen(input,"r")) ==NULL) {
    __OUT_BASE_LINE__("can not open altimetric data file: %s\n",input);
    exit(-1);
    }

  n=999;
  for(i=0;i<34;i++) fgets(line,n,in);
  nitems=fscanf(in,"%d",&nmes);
  fgets(line,n,in);
  fgets(line,n,in);

  printf("file %s : %d \n",input,nmes);

  metadata.data=new data_t[nmes];
  buffer=metadata.data;

  printf("pas de test interior car on suppose toutes les donnes valides\n");
  for(k=0;k<nmes;k++) {
      nitems=fscanf(in,"%d %lf %lf %lf",&cycle, &lon, &lat, &(buffer[k].time));
      if(nitems != 4) {__OUT_BASE_LINE__("read error exit \n");exit(2);}
      buffer[k].lon=lon;
      buffer[k].lat=lat;
      buffer[k].cycle=cycle;

      for(j=0;j<26;j++) {
          nitems=fscanf(in,"%f",&(buffer[k].values[j]));
          if(nitems!=1) {__OUT_BASE_LINE__("read error exit \n");exit(2);}
        }

//       if(npolygones !=0)
// 	{
// 	  if(lon>180) lon-=360.0;
// 	  interior=plg_TestInterior(lon,lat,polygones,npolygones);
// 	  if(interior !=1 ) continue;
// 	}

    }
  fclose(in);

  metadata.count=nmes;
  metadata.mask=99.9999;

  printf("%s starts %s, finishes %s \n",input,sgetcnesdate(buffer[0].time*24.),sgetcnesdate(buffer[nmes-1].time*24.));

/*
   IB = mog2D-G
   TIDES = FES2002
 */

  keep=(int *)calloc(nmes,sizeof(int));  //calloc necessary    not initialised
  serie=(double *)malloc(nmes*sizeof(double));
  t=(double *)malloc(nmes*sizeof(double));

/*-----------------------------------------------------------------------------
  data rejection when masked correction */

  for(k=0,n=0;k<nmes;k++) {
      if(buffer[k].values[0]==metadata.mask) {keep[k]=1;n++;continue;}
      if(buffer[k].values[MOG]==metadata.mask){keep[k]=1;n++;continue;}
      if(buffer[k].values[10]==metadata.mask) {keep[k]=1;n++;continue;}
      if(buffer[k].values[11]==metadata.mask) {keep[k]=1;n++;continue;}
      if(buffer[k].values[TIDE]==metadata.mask) {keep[k]=1;n++;continue;}
      if(buffer[k].values[24]==metadata.mask) {keep[k]=1;n++;continue;}
    }

  nmes=decale_metadata(buffer,keep,nmes);


  /* internal test */
  for(k=0;k<nmes;k++) {
      tides=buffer[k].values[10]+buffer[k].values[11]+buffer[k].values[TIDE];
      mss=buffer[k].values[24];
      serie[k]=buffer[k].values[0]-buffer[k].values[MOG]-tides-mss;
    }
  sdum=get_statistics(serie, metadata.mask, nmes);


/*-----------------------------------------------------------------------------
  3-sigma low-pass filter */

  z=(double *)malloc(nmes*sizeof(double));
  for(k=0;k<nmes;k++) keep[k]=0;
  threshold=3;

  for(k=0,n=0;k<nmes;k++) {
      z[k]=(serie[k]-sdum.mean);
      if(z[k]*z[k] > threshold*sdum.std*sdum.std ) {keep[k]=1;n++;}
    }

  nmes=decale_metadata(buffer,keep,nmes);

  /* internal test */
  for(k=0;k<nmes;k++) {
      tides=buffer[k].values[10]+buffer[k].values[11]+buffer[k].values[TIDE];
      mss=buffer[k].values[24];
      serie[k]=buffer[k].values[0]-buffer[k].values[MOG]-tides-mss;
    }
  sdum=get_statistics(serie, metadata.mask, nmes);

/*-----------------------------------------------------------------------------
  variables description
  0  is T/P SSH without masked value nor absurd values (>3sigma)
  28 is 0 low-passed filtered  (Loess)
  29 is 28 Sa filtered (harmonic analysis) */


 /*-----------------------------------------------------------------------------
    Loess filter each T/P pass*/

  for(k=0;k<nmes;k++) serie[k]=metadata.mask;
  i=0;
  j=0;
  do
    {
      cycle=buffer[npoints].cycle;
      while(buffer[npoints].cycle==cycle) npoints++;

      switch (npoints-i) {
        case (1):
          tides=buffer[i+k].values[10]+buffer[i+k].values[11]+buffer[i+k].values[TIDE];
          mss=buffer[i+k].values[24];
          buffer[i+k].values[28]=buffer[i+k].values[0]-buffer[i+k].values[MOG]-tides-mss;
          break;
        
        default:
            /* cut-off frequency: 30km -> 'cut' s */
            dum=geo_km(buffer[i].lon,buffer[i].lat,buffer[npoints-1].lon,buffer[npoints-1].lat);
            if(dum==0.0) {
                printf("distance error cycle %d\n",cycle);
                continue;
              }
            cut=30*(buffer[i].time-buffer[npoints-1].time)/dum;
            if (cut<0) cut*=-1;
          
            for(k=0;k<npoints-i;k++) {
                tides=buffer[i+k].values[10]+buffer[i+k].values[11]+buffer[i+k].values[TIDE];
                mss=buffer[i+k].values[24];
                z[k]=buffer[i+k].values[0]-buffer[i+k].values[MOG]-tides-mss;
              }
            for(k=0;k<npoints-i;k++) t[k]=buffer[i+k].time;
            status=loess1d_irregular(npoints-i,cut,z,t,(double) metadata.mask,serie);
            for(k=0;k<npoints-i;k++) buffer[i+k].values[28]=serie[k];
            for(k=0;k<npoints-i;k++) buffer[i+k].values[29]=serie[k];
        }/* end switch*/
      i=npoints;
    }
  while(npoints<nmes);

  /* internal test */
  for(k=0;k<nmes;k++)  serie[k]=buffer[k].values[28];
  sdum=get_statistics(serie, metadata.mask, nmes);

  for(k=0;k<nmes;k++) {
      tides=buffer[k].values[10]+buffer[k].values[11]+buffer[k].values[TIDE];
      mss=buffer[k].values[24];
      buffer[k].values[28]+=(buffer[k].values[MOG]+tides+mss);
      buffer[k].values[29]+=(buffer[k].values[MOG]+tides+mss);
    }


/*-----------------------------------------------------------------------------
  remove annual signal*/

  status=spectrum_init(&s);
  status=addwave(&s, wSa);
  spectrum_terminate(s);

  residuals=(double *)malloc(nmes*sizeof(double));

  for(k=0;k<nmes;k++) {
      t[k]=buffer[k].time;
      serie[k]=buffer[k].values[29];
    }
  status=harmonic_analysis_old(serie,t,residuals,&mean,nmes, s);
  for(k=0;k<nmes;k++) buffer[k].values[29]=mean+residuals[k]+buffer[k].values[24];
  sdum=get_statistics(residuals, (double) metadata.mask, nmes);
 
  metadata.data=buffer;
  metadata.count=nmes;

/*-----------------------------------------------------------------------------
   writing T/P measurements locations (.nod file) */

  sprintf(nodefile,"%s.nod",input);
  status= savenodes(nodefile,  metadata);

  
  free(serie); serie=NULL;
  free(t);t=NULL;
  free(residuals);residuals=NULL;
  free(z);z=NULL;
  free(keep);keep=NULL;

  return(metadata);

}


/*----------------------------------------------------------------------------*/

void write_list_header(FILE *out,int n)
     
/*----------------------------------------------------------------------------*/
{

  fprintf(out,"#-- HEADER -------------------------------------------\n");
  fprintf(out,"# Column 1 : date in days referred to\n");
  fprintf(out,"#            CNES date (01-JAN-1950 00:00:00.0)\n");
  fprintf(out,"# Column 2 : sea surface height (in meters) - Mog2DG -Loading -Solid ....\n");
  fprintf(out,"#-- HEADER END ---------------------------------------\n");
  fprintf(out,"Number of crossover points :\n");
  fprintf(out,"%d\n",n);
  fprintf(out,"#-----------------------------------------------------\n");

  fflush(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void save_meantracks_CLS(serie_t TPdata,char *track,float *mss,plg_t *polygones,int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,l,n,nitems;
  int nXpoints,dimD;
  int idum,nmesmax=1000;
  double x,y;
  double *t,*h,*residuals,mean;
  char *sdum,line[300],a;
  FILE *meantrk_file;
  serie_t *CLSprojected,*projected;
  spectrum_t s;
  statistic_t statdum;
  
 
  printf("Creating mean tracks time series...\n");
  dimD=TPdata.count;

  /*---------------------------------------------------------------------
    reading CLS mean tracks time series */

  sprintf(line,"../data/track-ref.CLS.TP.%s.dat",track);

  if( (meantrk_file=fopen(line,"r")) == NULL ) {
      sprintf(line,"can not open ../data/track-ref.TP.%s.dat",track);
      gmerror(line);
    }

  do  fgets(line,sizeof(line),meantrk_file);
  while (strncmp(line,"#-- HEADER END ---------------------------------------",20) != 0 );
  
  fgets(line,sizeof(line),meantrk_file);
  fgets(line,sizeof(line),meantrk_file);
  sscanf(line, "%d", &nXpoints);
    
  CLSprojected=(serie_t *)malloc(nXpoints*sizeof(serie_t));

  for(i=0;i<nXpoints;i++) {
      fgets(line,sizeof(line),meantrk_file);
      if( (nitems=fscanf(meantrk_file,"Pt  :%d\n",&idum)) != 1) gmerror("error in reading Xover number");
      if( (nitems=fscanf(meantrk_file,"lon :%lf\n",&x)) != 1) gmerror("error in reading Xover lon");
      if( (nitems=fscanf(meantrk_file,"lat :%lf\n",&y)) != 1) gmerror("error in reading Xover lat");
      if( (nitems=fscanf(meantrk_file,"Mes :%d",&(CLSprojected[i].count))) != 1) gmerror("error in reading Xover Mes");
      
      CLSprojected[i].mask=99.9999;
      CLSprojected[i].data=new data_t[CLSprojected[i].count];

      for(j=0;j<CLSprojected[i].count;j++) {
          CLSprojected[i].data[j].cycle=idum;
          CLSprojected[i].data[j].lon=x;
          CLSprojected[i].data[j].lat=y;

          nitems=fscanf(meantrk_file,"%lf",&(CLSprojected[i].data[j].time));
          for(k=0;k<24;k++) {
              nitems=fscanf(meantrk_file,"%f",&(CLSprojected[i].data[j].values[k]));
              if(nitems!=1)  gmerror("error in reading Xover values");
            }

          /*---------------------------------------------------------------------
            assuming geophysical and MSS corrections already made by CLS */

          CLSprojected[i].data[j].values[28]=
            CLSprojected[i].data[j].values[0]     /* SLA */
            -CLSprojected[i].data[j].values[1]    /* MOG2D-G */
            -CLSprojected[i].data[j].values[10]   /* loading effects, not included in CLS Xover */
            -CLSprojected[i].data[j].values[11]   /* solid Earth tide, not included in CLS Xover */
            -CLSprojected[i].data[j].values[12];  /* MOG2D-medsea (tides) */
        
          do {a=fgetc(meantrk_file);} while ( a!= '\n');
        }
    }/* end for i */
  fclose(meantrk_file);


 
  /*---------------------------------------------------------------------
    compute T/P ssh time series on CLS cross over locations
    method: substracting a mean sea surface to T/P along track data
             in a 0.025deg circle around each cross over location  */

  projected=(serie_t *)malloc(nXpoints*sizeof(serie_t));
  for(i=0;i<nXpoints;i++) {
      projected[i].data=new data_t[nmesmax];
      projected[i].count=0;
      projected[i].mask=99.9999;
    }

 for(i=0;i<nXpoints;i++) {
      k=0;
      for(j=0;j<dimD;j++) {
          x=TPdata.data[j].lon;
          y=TPdata.data[j].lat;
          if( (fabs(x-CLSprojected[i].data[0].lon) < 0.025) && (fabs(y-CLSprojected[i].data[0].lat) < 0.025) ) {
              projected[i].count++;
              projected[i].data[k].cycle=CLSprojected[i].data[0].cycle;
              projected[i].data[k].lon=CLSprojected[i].data[0].lon;
              projected[i].data[k].lat=CLSprojected[i].data[0].lat;
              projected[i].data[k].time=TPdata.data[j].time;

              projected[i].data[k].values[0]= TPdata.data[j].values[0]-mss[j];
              for(l=1;l<24;l++) projected[i].data[k].values[l]=TPdata.data[j].values[l];
/* 	      projected[i].data[k].values[28]= TPdata.data[j].values[28]-mss[j];     !!!!!  */
              projected[i].data[k].values[28]= TPdata.data[j].values[28]-TPdata.data[j].values[12]-mss[j];
              k++;
            }
        }
    }/* loop on XOver */

 sprintf(line,"track-ref.POC.TP.%s.dat",track);
 if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track file");

 i=0;
 for(n=0;n<nXpoints;n++) {
     if(projected[n].count==0) continue;
     
     x=CLSprojected[n].data[0].lon;
     y=CLSprojected[n].data[0].lat;
     idum=plg_TestInterior(x,y,polygones,npolygones);
     if(idum == 1) i++;
   }
 write_list_header(meantrk_file,i);

 idum=0;
 for(n=0;n<nXpoints;n++) {
     if(projected[n].count==0) continue;
     
     x=CLSprojected[n].data[0].lon;
     y=CLSprojected[n].data[0].lat;
     idum=plg_TestInterior(x,y,polygones,npolygones);
     if(idum == 1) {
 
         /*-----------------------------------------------------------------------------
           remove annual signal to POC mean tracks time series */

/* 	 idum=spectrum_init(&s); */
/* 	 idum=addwave(&s, wSa); */
/* 	 spectrum_terminate(s); */
 
/* 	 h=malloc(projected[n].count*sizeof(double)); */
/* 	 t=malloc(projected[n].count*sizeof(double)); */
/* 	 residuals=malloc(projected[n].count*sizeof(double)); */

/* 	 for(k=0;k<projected[n].count;k++) */
/* 	   { */
/* 	     t[k]=projected[n].data[k].time; */
/* 	     h[k]=(double) projected[n].data[k].values[28]; */
/* 	   } */

/* 	 status=harmonic_analysis(h,t,residuals,&mean,projected[n].count, s); */

/* 	 for(k=0;k<projected[n].count;k++) projected[n].data[k].values[29]=residuals[k]+mean; */
/* 	 for(k=0;k<projected[n].count;k++) h[k]=(double) projected[n].data[k].values[29]; */
/* 	 statdum=get_statistics(h,projected[n].mask,projected[n].count); */
 
         /*-----------------------------------------------------------------------------
           writing mean tracks time series */

         fprintf(meantrk_file,"Pt  : %d\n",projected[n].data[0].cycle);
         fprintf(meantrk_file,"lon : %lf\n",projected[n].data[0].lon);
         fprintf(meantrk_file,"lat : %lf\n",projected[n].data[0].lat);
         fprintf(meantrk_file,"Mes : %d\n",projected[n].count);
         
         for(k=0;k<projected[n].count;k++) {
             fprintf(meantrk_file,"%12.6f",projected[n].data[k].time);
             for(l=0;l<24;l++) fprintf(meantrk_file," %7.4f",projected[n].data[k].values[l]);
/* 	     fprintf(meantrk_file," %7.4f %7.4f\n",projected[n].data[k].values[28],projected[n].data[k].values[29]); */
             fprintf(meantrk_file,"\n");
           }
         fprintf(meantrk_file,"#-----------------------------------------------------\n");
 
/* 	 free(h); */
/* 	 free(t); */
/* 	 free(residuals); */
 
         /*-----------------------------------------------------------------------------
           remove annual signal to CLS mean tracks time series */
 
/* 	 h=malloc(CLSprojected[n].count*sizeof(double)); */
/* 	 t=malloc(CLSprojected[n].count*sizeof(double)); */
/* 	 residuals=malloc(CLSprojected[n].count*sizeof(double)); */
 
/* 	 for(k=0;k<CLSprojected[n].count;k++) h[k]=(double) CLSprojected[n].data[k].values[0]; */
/* 	 statdum=get_statistics(h,CLSprojected[n].mask,CLSprojected[n].count); */

/* 	 for(k=0;k<CLSprojected[n].count;k++) */
/* 	   { */
/* 	     t[k]=CLSprojected[n].data[k].time; */
/* 	     h[k]=(double) CLSprojected[n].data[k].values[28];	 */
/* 	   } */

/* 	 status=harmonic_analysis(h,t,residuals,&mean,CLSprojected[n].count, s); */

/* 	 for(k=0;k<CLSprojected[n].count;k++) CLSprojected[n].data[k].values[29]=residuals[k]+mean; */
/* 	 for(k=0;k<CLSprojected[n].count;k++) h[k]=(double) CLSprojected[n].data[k].values[29]; */
/* 	 statdum=get_statistics(h,CLSprojected[n].mask,CLSprojected[n].count); */
 
         /*-----------------------------------------------------------------------------
           writing CLS mean tracks time series */
/* 	 if( (meantrk_file=fopen("out2","w")) == NULL ) gmerror("can not open /calcul/maewo/roblou/out2"); */

/* 	 write_list_header(meantrk_file); */
/* 	 fprintf(meantrk_file,"Pt  : %d\n",CLSprojected[n].data[0].cycle); */
/* 	 fprintf(meantrk_file,"lon : %lf\n",CLSprojected[n].data[0].lon); */
/* 	 fprintf(meantrk_file,"lat : %lf\n",CLSprojected[n].data[0].lat); */
/* 	 fprintf(meantrk_file,"Mes : %d\n",CLSprojected[n].count); */

/* 	 for(k=0;k<CLSprojected[n].count;k++) */
/* 	   { */
/* 	     fprintf(meantrk_file,"%12.6f",CLSprojected[n].data[k].time);  */
/* 	     for(l=0;l<24;l++) fprintf(meantrk_file," %7.4f",CLSprojected[n].data[k].values[l]);  */
/* 	     fprintf(meantrk_file,"\n");  */
/* 	     fprintf(meantrk_file," %7.4f %7.4f\n",CLSprojected[n].data[k].values[28],CLSprojected[n].data[k].values[29]);  */
/* 	   } */
/* 	 fprintf(meantrk_file,"#-----------------------------------------------------\n"); */
/* 	 free(h); */
/* 	 free(t); */
/* 	 free(residuals); */
     
       } /* if idum */

   }   /* loop on nXpoints */

 fclose(meantrk_file);

} /* end */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


void save_meantracks_CTOH(serie_t TPdata,char *alti_file, int pair,float *mss,grid_t *nominal_track)


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int i,j,k,l,n,nitems,one_time=1;
  int nXpoints,dimD,*keep,OK;
  int idum,nmesmax=1000;
  float size;
  double x,y,*lon,*lat,ref,d[2];
  char *sdum,line[300];
  FILE *meantrk_file;
  serie_t *projected;
  statistic_t stat;
  int nb_cycles,cycle_min=9999;


  printf("Creating mean tracks time series...\n");
  dimD=TPdata.count;


  nXpoints=nominal_track->nx;
  keep=(int *)calloc(nXpoints,sizeof(int));
    
  for(i=0;i<nXpoints;i++) if(nominal_track->x[i]>180) nominal_track->x[i]-=360;


  /*---------------------------------------------------------------------
    organizing locations to be along track */

  if(pair%2==0) /* descending T/P, J1, or T/P interleaved track!, must flip indexes...*/
    for(i=0,j=nXpoints-1;i<(int) floor(nXpoints/2.0);i++,j--) {
        x=nominal_track->x[i];
        nominal_track->x[i]=nominal_track->x[j];
        nominal_track->x[j]=x;

        y=nominal_track->y[i];
        nominal_track->y[i]=nominal_track->y[j];
        nominal_track->y[j]=y;
      }

  /*---------------------------------------------------------------------
    compute T/P ssh time series on CTOH mean tracks locations
    method: substracting a mean sea surface to T/P along track data */
  for(i=0;i<dimD-1;i++)
    TPdata.data[i].used=0;
  
  nb_cycles=1;
  for(i=0;i<dimD-1;i++) {
      if( TPdata.data[i+1].cycle-TPdata.data[i].cycle>=1) nb_cycles++;
      cycle_min=MIN(cycle_min,TPdata.data[i].cycle);
    }
  printf("%d cycles\n",nb_cycles);
  printf("premier cycle : %d\n",cycle_min);

  //methode de thierry .....

   projected=(serie_t *)calloc(nXpoints,sizeof(serie_t));   //attention calloc obligatoire car pas d'initialisation !!!!!
   for(i=0;i<nXpoints ;i++)projected[i].data=(data_t *)calloc(nb_cycles,sizeof(data_t));
   //projection sur la trace nominale ...

   //Je calcule le vecteur entre le point 1 et le point 2 de la trace nominale ...
   typedef struct
    {
      double  x,y;
      double norme;
    } vector_t;
   vector_t dir,dat;
   double prod_scal;
  for(i=0;i<nXpoints-1 ;i++)//attention probleme limite a la fin de la serie ....
    {
      dir.x=(nominal_track->x[i+1]-nominal_track->x[i]);  //attention probleme limite a la fin de la serie ....
      dir.y=(nominal_track->y[i+1]-nominal_track->y[i]);
      dir.norme=sqrt(dir.x*dir.x+dir.y*dir.y);
      //normalisation
      dir.x/=dir.norme;
      dir.y/=dir.norme;
    
      // je calcule le vecteur entre le point 1 nominale et les points de mesure du cycle
      size=0.5; //taille d'une boite pour optimiser la recherche
      for(j=0;j<dimD;j++)
        if( (fabs(TPdata.data[j].lon-nominal_track->x[i]) < size) && (fabs(TPdata.data[j].lat-nominal_track->y[i]) < size) ) {	
           //dat.x=(nominal_track->x[i]-TPdata.data[j].lon);
           //dat.y=(nominal_track->y[i]-TPdata.data[j].lat);
           dat.x= geo_km(nominal_track->x[i],nominal_track->y[i],TPdata.data[j].lon,nominal_track->y[i]);
           dat.y= geo_km(nominal_track->x[i],nominal_track->y[i],nominal_track->x[i],TPdata.data[j].lat);
           dat.norme=sqrt(dat.x*dat.x+dat.y*dat.y);
           //calcul du produit scalaire
           prod_scal=dir.x*dat.x+dir.y*dat.y;
           d[0]=geo_km(nominal_track->x[i],nominal_track->y[i],nominal_track->x[i+1],nominal_track->y[i+1]);
           //attribution ou non de la donnee au point nominal ...
           if( (fabs(prod_scal)<=0.5*d[0])&&(TPdata.data[j].used!=1) )// --> on prend la donnee si elle n'est pas deja prise ailleurs !!!;
             {
               idum=projected[i].count;
               projected[i].data[idum].cycle=TPdata.data[j].cycle;
               projected[i].data[idum].lon  =TPdata.data[j].lon;
               projected[i].data[idum].lat  =TPdata.data[j].lat;
               projected[i].data[idum].time =TPdata.data[j].time;
               projected[i].data[idum].values[0]= TPdata.data[j].values[0]-mss[j];
               for(l=1;l<24;l++) projected[i].data[idum].values[l]=TPdata.data[j].values[l];
               projected[i].data[idum].values[28]= TPdata.data[j].values[28]-mss[j];
               projected[i].data[idum].values[29]= TPdata.data[j].values[29]-mss[j];
               TPdata.data[j].used=1;
               projected[i].count++;
             }
   
         }
    }




  //methode de laurent .....

    
  /*---------------------------------------------------------------------
    1st step: collecting data in square boxes */

//   size=0.026;	
//   projected=(serie_t *)malloc(nXpoints*sizeof(serie_t));
//   for(i=0;i<nXpoints;i++)
//     {
//       projected[i].data=new data_t[nmesmax];
//       projected[i].count=0;
//       projected[i].mask=99.9999;
//     }
//   keep=(int *)realloc(keep,dimD*sizeof(int));
//   for(j=0;j<dimD;j++) keep[j]=999;

//   for(i=0;i<nXpoints;i++)
//     {
//       idum=-1;
//       k=0;
//       for(j=0;j<dimD;j++)
//      {
//        x=TPdata.data[j].lon;
//        y=TPdata.data[j].lat;
//        if( (fabs(x-nominal_track->x[i]) < size) && (fabs(y-nominal_track->y[i]) < size) )
// 	 {
// 	   projected[i].count++;
// 	   projected[i].data[k].cycle=TPdata.data[j].cycle;
// 	   projected[i].data[k].lon=x;
// 	   projected[i].data[k].lat=y;

// 	   projected[i].data[k].time=TPdata.data[j].time;
// 	   projected[i].data[k].values[0]= TPdata.data[j].values[0]-mss[j];	
// 	   for(l=1;l<24;l++) projected[i].data[k].values[l]=TPdata.data[j].values[l];
// 	   projected[i].data[k].values[28]= TPdata.data[j].values[28]-mss[j];
// 	   projected[i].data[k].values[29]= TPdata.data[j].values[29]-mss[j];
// 	   keep[j]=1;

// 	   k++;
// 	 }
//      }
//     }/* loop on raw data */



//   idum=0;
//   for(i=0;i<dimD;i++) if(keep[i]==999) idum++;
//   if(idum==0) printf("all along track data selected\n");
//   else printf("%d along track data not in selection\n",idum);

//   /*---------------------------------------------------------------------
//     2nd step: assigning nearest data each cycle to each ref location*/

//   for(i=0;i<nXpoints;i++)
//     {
//       idum=projected[i].data[0].cycle;
//       n=0;
//       for(k=0;k<projected[i].count;)
//      {
//        OK=1;
//        one_time=1;
//        ref=999;
//        j=k;
//        do{j++;} while(projected[i].data[j].cycle==idum) ;

//        for(nitems=k;nitems<j;nitems++)
// 	 {
// 	   d[0]=geo_km(nominal_track->x[i],nominal_track->y[i],projected[i].data[nitems].lon,projected[i].data[nitems].lat);
// 	   if(i<nXpoints)
// 	     d[1]=geo_km(nominal_track->x[i+1],nominal_track->y[i+1],projected[i].data[nitems].lon,projected[i].data[nitems].lat);
// 	   else d[1]=999;

// 	   ref=MIN(ref,d[0]);
// 	   ref=MIN(ref,d[1]);
           
// 	   l=0;
// 	   if(i>0) /* test if the data is already choosen */
// 	     while(projected[i-1].data[l].time<=projected[i].data[nitems].time)
// 	       {
// 		 if(l>projected[i-1].count) break;
// 		 if(projected[i].data[nitems].time==projected[i-1].data[l].time)  OK=0;
// 		 l++;
// 	       }

// 	   if( (ref==d[0]) && (OK) )
// 	     {
// 	       projected[i].data[n].cycle=projected[i].data[nitems].cycle;
// 	       projected[i].data[n].time=projected[i].data[nitems].time;
// 	       projected[i].data[n].values[0]=projected[i].data[nitems].values[0];	
// 	       for(l=1;l<24;l++) projected[i].data[n].values[l]=projected[i].data[nitems].values[l];
// 	       projected[i].data[n].values[28]=projected[i].data[nitems].values[28];
// 	       projected[i].data[n].values[29]=projected[i].data[nitems].values[29];
// 	       one_time=0;
// 	     }

// 	   OK=1;
// 	 }
//        idum=projected[i].data[j].cycle;
//        k=j;
//        if(one_time==0) n++;
//      }
//       projected[i].count=n;
//     }/* loop on ref time series */

  idum=0;
  for(i=0;i<nXpoints;i++) idum+=projected[i].count;
  printf("total data projected: %d (%d)\n",idum,dimD);

  idum=0;
  for(i=0;i<nXpoints;i++)
    for(k=0;k<projected[i].count-1;k++) {
        if( projected[i].data[k+1].time-projected[i].data[k].time<1)
          printf("pb serie %d k=%d (count:%d)\n",i+1,k+1,++idum);
      }

  idum=0;
  for(n=0;n<nXpoints;n++) if(projected[n].count==0) idum++;

  /*-----------------------------------------------------------------------------
    writing mean tracks time series (track-ref.POC.TP.*.dat file) */

  sprintf(line,"%s.Sa_not_removed_analyse_list",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      i++;
      if(projected[n].count==0) continue;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",nominal_track->x[n]);
      fprintf(meantrk_file,"lat : %lf\n",nominal_track->y[n]);
      fprintf(meantrk_file,"Mes : %d\n",projected[n].count);
         
      for(k=0;k<projected[n].count;k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k].time);
          fprintf(meantrk_file," %7.4f",projected[n].data[k].values[28]-projected[n].data[k].values[1]-projected[n].data[k].values[10]-projected[n].data[k].values[11]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);

  /*-----------------------------------------------------------------------------

  /*-----------------------------------------------------------------------------
    mapping mean tracks time series location ( *.nod file)*/

  sprintf(line,"track-ref.POC.%s.nod",alti_file);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean location file");
  //fprintf(meantrk_file,"%d\n0\n%d\n",nXpoints-idum,nXpoints-idum);
  fprintf(meantrk_file,"%d\n0\n%d\n",nXpoints,nXpoints);
  for(n=0;n<nXpoints;n++) {
      //if(projected[n].count==0) continue;
      fprintf(meantrk_file,"%7.4f %7.4f %d\n",nominal_track->x[n],nominal_track->y[n],projected[n].count);
    }
  fclose(meantrk_file);meantrk_file=NULL;

 
  free(nominal_track->x);nominal_track->x=NULL;
  free(nominal_track->y);nominal_track->y=NULL;
  free(keep);keep=NULL;
  for(i=0;i<nXpoints;i++){  free(projected[i].data); projected[i].data=NULL;}
  free(projected); projected=NULL;


} /* end */

