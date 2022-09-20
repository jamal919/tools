 #include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "netcdf-proto.h"
#include "poc-netcdf.def"
#  define RANK_mss 2
#  define RANK_cls 2
#  define RANK_innovation 2
#  define RANK_geoid 2
#  define RANK_topo 2
  /* IB = mog2D medsea */
#define MOG 4
#define TIDE 12 /* regional model */
#include "poc-time.h"
#include "tides.h"
#include "tides.def"
#include "filter.h"
#include "statistic.h"
#include "polygones.h"
#include "geo.h"
#include "map.h"
#include "rutin.h"


   /* variable ids */
   int lon_id;
   int lat_id;
   int mss_id;
   int cls_id;
   int innovation_id;
   int geoid_id;
   int topo_id;



//theses functions are defined is this CPP file
extern int mss_createfile(char *filename,size_t x_len, size_t y_len, grid_t grid);
extern int writefile(char *filename, float *mss, float *cls, float *ino);
extern int savenodes(char *filename, serie_t metadata);
extern serie_t load_metadata_ref(char *input, int point);
extern int decale_metadata(data_t *data,int *keep,int nmes);
extern serie_t load_metadata_raw(char *input,  plg_t *polygones, int npolygones, char *zone);
extern void write_list_header(FILE *out,int n);
extern void save_meantracks_CTOH(serie_t TPdata,char *sat, char *track,float *mss,plg_t *polygones,int npolygones,char *root);

/*------------------------------------------------------------------------
  
  Nom         :  createfile

  --------------------------------------------------------------------------*/
/* *------------------------------------------------------------------------
  attribute name "axis" changed in "content" to comply with CF standard*/
int mss_createfile(char *filename,size_t x_len, size_t y_len, grid_t grid)

{		
   int  ncid;			/* netCDF id */

   /* dimension ids */
   int x_dim;
   int y_dim;

   /* variable ids
   int lon_id;
   int lat_id;
   int mss_id;
   int cls_id;
   int innovation_id;
   int geoid_id;
   int topo_id;*/


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
  spectrum_t s,solved;
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
      if(buffer[k].values[0] != metadata.mask) {
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
  double mask=1.e+35;
  harmonic_analysis_with_parsing(serie,mask,t,residuals,&mean,count, s,solved, 0);
 
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

  free(serie);
  free(t);
  free(residuals);

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

serie_t load_metadata_raw(char *input,  plg_t *polygones, int npolygones, char *zone)

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
  spectrum_t s,solved;
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

  k=0;
  while (!feof(in)) {
    nitems=fscanf(in,"%d %lf %lf %lf",&cycle, &lon, &lat, &(buffer[k].time));
    if(nitems != 4) break;
    buffer[k].lon=lon;
    buffer[k].lat=lat;
    buffer[k].cycle=cycle;

    for(j=0;j<26;j++) {
      nitems=fscanf(in,"%f",&(buffer[k].values[j]));
      if(nitems!=1) break;
      }

    if(npolygones !=0) {
      interior=plg_TestInterior(lon,lat,polygones,npolygones);
      if(interior !=1 ) continue;
      }

    k++;
    }
  nmes=k;
  fclose(in);

  metadata.count=nmes;
  metadata.mask=99.9999;

  printf("%s starts %s, finishes %s \n",input,sgetcnesdate(buffer[0].time*24.),sgetcnesdate(buffer[nmes-1].time*24.));

/*
   IB = mog2D-G
   TIDES = FES2002
 */

  keep=(int *)calloc(nmes,sizeof(int));
  serie=(double *)malloc(nmes*sizeof(double));
  t=(double *)malloc(nmes*sizeof(double));

/*-----------------------------------------------------------------------------
  data rejection when masked correction */

  for(k=0;k<nmes;k++) keep[k]=0;
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
  do {
      cycle=buffer[npoints].cycle;
      while(buffer[npoints].cycle==cycle) npoints++;

      switch (npoints-i) {
        case 1:
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
          break;
        }
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
  double mask=1.e+35;
  harmonic_analysis_with_parsing(serie,mask,t,residuals,&mean,nmes, s,solved, 0);
  for(k=0;k<nmes;k++) buffer[k].values[29]=mean+residuals[k]+buffer[k].values[24];

  sdum=get_statistics(residuals, (double) metadata.mask, nmes);
 
  metadata.data=buffer;
  metadata.count=nmes;

/*-----------------------------------------------------------------------------
   writing T/P measurements locations (.nod file) */

  sprintf(nodefile,"%s.%s.nod",input,zone);
  status= savenodes(nodefile,  metadata);

  
  free(serie);
  free(t);
  free(residuals);
  free(z);
  free(keep);

  return(metadata);

}


/*----------------------------------------------------------------------------*/

void write_list_header(FILE *out,int n)
     
/*----------------------------------------------------------------------------*/
{

  fprintf(out,"#-- HEADER -------------------------------------------\n");
  fprintf(out,"# Column 1 : date in days referred to\n");
  fprintf(out,"#            CNES date (01-JAN-1950 00:00:00.0)\n");
  fprintf(out,"# Column 2 : sea surface height (in meters)\n");
  fprintf(out,"# Column 3 : MOG2D-G model sea level (in meters)\n");
  fprintf(out,"# Column 4 : MOG2D-G model S1 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 5 : MOG2D-G model S2 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 6 : MOG2D-medsea model sea level (in meters)\n");
  fprintf(out,"# Column 7 : MOG2D-medsea model S1 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 8 : MOG2D-medsea model S2 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 9 : inverted barometer (in meters)\n");
  fprintf(out,"# Column 10: inverted barometer S1 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 11: inverted barometer S2 atmospheric tide (in meters)\n");
  fprintf(out,"# Column 12: loading effects (in meters)\n");
  fprintf(out,"# Column 13: solid Earth tide (in meters)\n");
  fprintf(out,"# Column 14: MOG2D-medsea model ocean tide (in meters)\n");
  fprintf(out,"# Column 15: MOG2D-medsea model S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 16: MOG2D-medsea model S2 ocean tide (in meters)\n");
  fprintf(out,"# Column 17: FES99 model ocean tide (in meters)\n");
  fprintf(out,"# Column 18: FES99 model S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 19: FES99 model S2 ocean tide (in meters)\n");
  fprintf(out,"# Column 20: FES2002 model ocean tide (in meters)\n");
  fprintf(out,"# Column 21: FES2002 model S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 22: FES2002 model S2 ocean tide (in meters)\n");
  fprintf(out,"# Column 23: harmonic analysis ocean tide (in meters)\n");
  fprintf(out,"# Column 24: harmonic analysis S1 ocean tide (in meters)\n");
  fprintf(out,"# Column 25: harmonic analysis S2 ocean tide (in meters)\n");
  fprintf(out,"#-- HEADER END ---------------------------------------\n");
  fprintf(out,"Number of crossover points :\n");
  fprintf(out,"%d\n",n);
  fprintf(out,"#-----------------------------------------------------\n");

  fflush(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void save_meantracks_CTOH(serie_t TPdata,char *sat, char *track,float *mss,plg_t *polygones,int npolygones,char *root)

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
 
 
  /*---------------------------------------------------------------------
    temporary, no mean positions available for GFO */
  if(strcasecmp(sat,"GFO")==0) return;
  if(strcasecmp(sat,"ENV")==0) return;


  printf("Creating mean tracks time series...\n");
  dimD=TPdata.count;

  /*---------------------------------------------------------------------
    reading CTOH mean tracks positions */

  if( (strcasecmp(sat,"TP")==0)||(strcasecmp(sat,"J1")==0) ) {
      nXpoints=3127;

      lon=(double *) calloc(nXpoints,sizeof(double));
      lat=(double *) calloc(nXpoints,sizeof(double));
      keep=(int *)calloc(nXpoints,sizeof(int));

      sprintf(line,"/data/usrcto/produits/base/TOPEX/TRACES_NOMINALES/INDEX/lat.dat",track);
      if( (meantrk_file=fopen(line,"r")) == NULL ) gmerror("can not open CTOH mean track latitude file");
      for(i=0;i<nXpoints;i++) fscanf(meantrk_file,"%lf\n",&(lat[i]));
      fclose(meantrk_file);

      sprintf(line,"/data/usrcto/produits/base/TOPEX/TRACES_NOMINALES/INDEX/lon_%s.dat",track);
      if( (meantrk_file=fopen(line,"r")) == NULL ) gmerror("can not open  CTOH mean track longitude file");
      for(i=0;i<nXpoints;i++) fscanf(meantrk_file,"%lf\n",&(lon[i]));
      fclose(meantrk_file);
    }

  if(strcasecmp(sat,"TP2")==0) {
      nXpoints=3127;

      lon=(double *) calloc(nXpoints,sizeof(double));
      lat=(double *) calloc(nXpoints,sizeof(double));
      keep=(int *)calloc(nXpoints,sizeof(int));
  
      sprintf(line,"/data/usrcto/travail_en_cours/TRACES_NOMINALES_TOPEX_NO/lat.dat",track);
      if( (meantrk_file=fopen(line,"r")) == NULL ) gmerror("can not open CTOH mean track latitude file");
      for(i=0;i<nXpoints;i++) fscanf(meantrk_file,"%lf\n",&(lat[i]));
      fclose(meantrk_file);

      sprintf(line,"/data/usrcto/travail_en_cours/TRACES_NOMINALES_TOPEX_NO/lon_%s.dat",track);
      if( (meantrk_file=fopen(line,"r")) == NULL ) gmerror("can not open CTOH mean track longitude file");
      for(i=0;i<nXpoints;i++) fscanf(meantrk_file,"%lf\n",&(lon[i]));
      fclose(meantrk_file);
    }

 
  /*---------------------------------------------------------------------
    1/ rejecting positions out of local polygone */

  n=0;
  for(i=0;i<nXpoints;i++) {
      lon[i]=geo_recale(lon[i],0.,180.0);
  
      idum=plg_TestInterior(lon[i],lat[i],polygones,npolygones);
      if(idum!=1) {keep[i]=1;n++;}
    }

  for(i=0,j=0;i<nXpoints;i++)
    if(keep[i]!=1) {
        lon[j]=lon[i];
        lat[j]=lat[i];
        j++;
      }
  nXpoints-=n;

  for(i=0;i<nXpoints;i++) if(lon[i]>180) lon[i]-=360;


  /*---------------------------------------------------------------------
    organizing locations to be along track */

  j=atoi(track);
  if(j%2==0) /* descending T/P, J1, or T/P interleaved track!, must flip indexes...*/
    for(i=0,j=nXpoints-1;i<(int) floor(nXpoints/2.0);i++,j--) {
        x=lon[i];
        lon[i]=lon[j];
        lon[j]=x;

        y=lat[i];
        lat[i]=lat[j];
        lat[j]=y;
      }

  /*---------------------------------------------------------------------
    compute T/P ssh time series on CTOH mean tracks locations
    method: substracting a mean sea surface to T/P along track data */
  
  idum=1;
  for(i=0;i<dimD-1;i++)
    if( TPdata.data[i+1].cycle-TPdata.data[i].cycle>=1) idum++;
  printf("%d cycles\n",idum);
  nmesmax=idum*3;
      
  /*---------------------------------------------------------------------
    1st step: collecting data in square boxes */

  size=0.06;
  projected=(serie_t *)malloc(nXpoints*sizeof(serie_t));
  for(i=0;i<nXpoints;i++) {
      projected[i].data=new data_t[nmesmax];
      projected[i].count=0;
      projected[i].mask=99.9999;
    }
  keep=(int *)realloc(keep,dimD*sizeof(int));
  for(j=0;j<dimD;j++) keep[j]=999;

  for(i=0;i<nXpoints;i++) {
      idum=-1;
      k=0;
      for(j=0;j<dimD;j++) {
          x=TPdata.data[j].lon;
          y=TPdata.data[j].lat;
          if( (fabs(x-lon[i]) < size) && (fabs(y-lat[i]) < size) ) {	
              projected[i].count++;
              projected[i].data[k].cycle=TPdata.data[j].cycle;
              projected[i].data[k].lon=x;
              projected[i].data[k].lat=y;

              projected[i].data[k].time=TPdata.data[j].time;
              projected[i].data[k].values[0]= TPdata.data[j].values[0]-mss[j];
              for(l=1;l<24;l++) projected[i].data[k].values[l]=TPdata.data[j].values[l];
              projected[i].data[k].values[28]= TPdata.data[j].values[28]-mss[j];
              projected[i].data[k].values[29]= TPdata.data[j].values[29]-mss[j];
              keep[j]=1;

              k++;
            }
        }
    }/* loop on raw data */

  idum=0;
  for(i=0;i<dimD;i++) if(keep[i]==999) idum++;
  if(idum==0) printf("all along track data selected\n");
  else printf("%d along track data not in selection\n",idum);

  /*---------------------------------------------------------------------
    2nd step: assigning nearest data each cycle to each ref location*/

  for(i=0;i<nXpoints;i++) {
      idum=projected[i].data[0].cycle;
      n=0;
      for(k=0;k<projected[i].count;) {
          OK=1;
          one_time=1;
          ref=999;
          j=k;
          do{j++;} while(projected[i].data[j].cycle==idum) ;

          for(nitems=k;nitems<j;nitems++) {
              d[0]=geo_km(lon[i],lat[i],projected[i].data[nitems].lon,projected[i].data[nitems].lat);
              if(i<nXpoints)
                d[1]=geo_km(lon[i+1],lat[i+1],projected[i].data[nitems].lon,projected[i].data[nitems].lat);
              else d[1]=999;

              ref=MIN(ref,d[0]);
              ref=MIN(ref,d[1]);
              
              l=0;
              if(i>0) /* test if the data is already choosen */
                while(projected[i-1].data[l].time<=projected[i].data[nitems].time)
                  {
                    if(l>projected[i-1].count) break;
                    if(projected[i].data[nitems].time==projected[i-1].data[l].time)  OK=0;
                    l++;
                  }

              if( (ref==d[0]) && (OK) ) {
                  projected[i].data[n].cycle=projected[i].data[nitems].cycle;
                  projected[i].data[n].time=projected[i].data[nitems].time;
                  projected[i].data[n].values[0]=projected[i].data[nitems].values[0];
                  for(l=1;l<24;l++) projected[i].data[n].values[l]=projected[i].data[nitems].values[l];
                  projected[i].data[n].values[28]=projected[i].data[nitems].values[28];
                  projected[i].data[n].values[29]=projected[i].data[nitems].values[29];
                  one_time=0;
                }

              OK=1;
            }
          idum=projected[i].data[j].cycle;
          k=j;
          if(one_time==0) n++;
        }
      projected[i].count=n;
    }/* loop on ref time series */

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

  sprintf(line,"track-ref.POC.%s.%s.dat",sat,track);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean track data file");
 
  write_list_header(meantrk_file,nXpoints-idum);
  for(n=0,i=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      i++;
      fprintf(meantrk_file,"Pt  : %d\n",i);
      fprintf(meantrk_file,"lon : %lf\n",lon[n]);
      fprintf(meantrk_file,"lat : %lf\n",lat[n]);
      fprintf(meantrk_file,"Mes : %d\n",projected[n].count);
         
      for(k=0;k<projected[n].count;k++) {
          fprintf(meantrk_file,"%12.6f",projected[n].data[k].time);
          fprintf(meantrk_file," %7.4f",projected[n].data[k].values[28]);
          for(l=1;l<24;l++) fprintf(meantrk_file," %7.4f",projected[n].data[k].values[l]);
          fprintf(meantrk_file,"\n");
        }
      fprintf(meantrk_file,"#-----------------------------------------------------\n");
    }   /* loop on nXpoints */
  fflush(meantrk_file);
  fclose(meantrk_file);

  /*-----------------------------------------------------------------------------
    mapping mean tracks time series location ( *.nod file)*/

  sprintf(line,"%s/track-ref.POC.%s.%s.nod",root,sat,track);
  if( (meantrk_file=fopen(line,"w")) == NULL ) gmerror("can not write in POC mean location file");
  fprintf(meantrk_file,"%d\n0\n%d\n",nXpoints-idum,nXpoints-idum);
  for(n=0;n<nXpoints;n++) {
      if(projected[n].count==0) continue;
      fprintf(meantrk_file,"%7.4f %7.4f %d\n",lon[n],lat[n],projected[n].count);
    }

  free(projected);
  for(i=0;i<nXpoints;i++) free(projected[i].data);
  free(lon);
  free(lat);
  free(keep);

} /* end */

