
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Sara Fleury        LEGOS/CNRS, Toulouse, France
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief Convert TUGOm output to simuSWOT-compliant landscape files.

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION

#include "config.h"

#include <stdio.h>
#include <unistd.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "tides.h"
#include "fe.h"
#include "map.h"
#include "topo.h"
#include "poc-time-classes.h"
#include "poc-netcdf-data.hpp"
#include "filter.h"
#include "statistic.h"
#include "ascii.h"
#include "archive.def"
#include "archive.h"
#include "functions.h"
#include "datastream.h"


#define ARCHIVE_FORMAT_ASCII  0
#define ARCHIVE_FORMAT_NETCDF 1
#define ARCHIVE_FORMAT_BINARY 2

#define GLOBAL      0
#define PER_YEAR    1
#define PER_MONTH   2
#define PER_WEEK    3
#define PER_DAY     4
#define PER_HOUR    5
#define PER_MINUTE  6
#define PER_SECOND  7
#define CUSTOM      8

extern int identify_BINARY(const char *file, const char *name, int *id, int verbose);
extern int stream_UGxBINARY(UGfield_t<double> & field, double t);

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int identify_BINARY(const char *file, const char *name, int *id, int verbose)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   if(strcmp(name,"elevation")==0) {
//     *id=0;
//     }
//   if(strcmp(name,"ubar")==0) {
//     *id=1;
//     }
//   if(strcmp(name,"vbar")==0) {
//     *id=2;
//     }
//   return(0);
// }


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int obc_ReadBinary(const char *datafile, double time, float **buffer, double *actual_time, double *next)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int step, status, d1, d2;
//   double start,shift;
// 
//   FILE *in;
//   date_t actual;
//   meta_archive_t info;
// 
// /**----------------------------------------------------------------------
//   to be fixed */
// //  getmodeldate(time, &actual);
// /**----------------------------------------------------------------------
//   to be fixed */
// //  printf("#obc_ReadBinary at %s in %s\n", sgetmodeldate(time), datafile);
// 
//   status = clx_archiveinfo(datafile, &info);
//   if(status != 0) {
//     printf("error while reading header in datafile: %s\n",datafile);
//     goto error;
//     }
// 
// /**---------------------------------------------------------------------
//   model time reference */
//   d1 = julian_day(t_reference.month, t_reference.day, t_reference.year);
//   
// /**---------------------------------------------------------------------
//   external archive time reference */
//   d2 = julian_day(info.reference.month, info.reference.day,info.reference.year);
// 
// /**---------------------------------------------------------------------
//   shift between model and archive time: archive time=model time -start */
//   shift = (d2 - d1) * 24 * 3600 - t_reference.second;
//   start = info.start + shift;
// 
//   step = 1 + (int) ((time - start) / info.sampling);
//   in = fopen(datafile, "r");
//   if(in == NULL) {
//     printf("error while opening the datafile %s\n", datafile);
//     goto error;
//     }
// 
//   status = clx_archiveread(in, info, step, buffer, actual_time);
//   if(status != 0) {
//     printf("error while reading the datafile at step =%d %s\n", step,datafile);
//     if(step < 1)
// /**----------------------------------------------------------------------
//   to be fixed */
// //      printf("model start before archive content : %s\n",sgetmodeldate(start));
//     goto error;
//     }
// 
//   fclose(in);
// 
// /**---------------------------------------------------------------------
//   convert archive time in model time */
//   *actual_time += shift;
//   *next = *actual_time +info.sampling;
//      
//   info.destroy();
// 
//   return (status);
// 
// error:
//   status = -1;
//   return (status);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int stream_UGxBINARY(UGfield_t<double> & field, double t)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// //   using namespace Keywords ;
// 
//   int status, k;
//   double time,next;
//   float *buffer[3];
//   filestream_t<float> *stream;
//   date_t actual, reference;
// 
//   stream=(filestream_t<float> *) field.stream;
//   if(stream==0) return(-1);
//   
//   stream->check(t);
//   
// /**----------------------------------------------------------------------
//   to be fixed */
// //  if(external_obc == 0)
// //    return (0);
// 
// /*------------------------------------------------------------------------------
//   scale factor for backward compatibility*/
// /**----------------------------------------------------------------------
//   to be fixed */
// //  if (tugo_cfg->boundaries->ArchivedBCUnits.actual_value == KEY_MKS) {
// //     obc_scale=1.0;
// //     }
// //   if (tugo_cfg->boundaries->ArchivedBCUnits.actual_value == KEY_CGS) {
// //     obc_scale=1.0e-02;
// //     }
// 
//   for(k = 0; k < 3; k++)
//     buffer[k] = new float[field.descriptor->nnodes];
//   
//   status = obc_ReadBinary(stream->filename.c_str(), t, buffer, &time, &next);
//   if(status!=0) return (-1);
//   
//   field.time=time;
//   field.next=next;
// /**----------------------------------------------------------------------
//   to be fixed */
// //  for(int n=0;n<field.descriptor->nnodes;n++) field.x[n]=obc_scale*buffer[stream->id][n];
// 
//   field.units=strdup("UNDOCUMENTED");
// 
//   for(k = 0; k < 3; k++)
//     delete[] buffer[k];
//   
//   return (0);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_copy(UGfield_t<double> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
{
  int n, status;
  
/*------------------------------------------------------------------------------
  tricky bypass, see below */
  sequence_t<UGfield_t<double> > *stream=(sequence_t<UGfield_t<double> > *)field.stream;
  
  status=stream->check(t);
  
  #pragma omp parallel for
  for(n=0;n<field.descriptor->nnodes;n++) {
    field.x[n]=stream->interpolate(t,n);
    }
  
  field.time=t;
//   field.next=stream->data.next;
  
  field.units =strdup(stream->frames[0].units);
  field.mask  =stream->frames[0].mask;
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_SGxUG(SGfield_t<double> & Sfield, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  
  interpolate an unstructured field on a structured grid

------------------------------------------------------------------------------*/
{
  int status;

  UGfield_t<double>  *stream;
  stream=(UGfield_t<double> *) Sfield.stream;
  
/*------------------------------------------------------------------------------
  update stream, i.e. unstructured field */
  status=stream->check(t);/* check because it may be used to interpolate bathymetry */
  
/*------------------------------------------------------------------------------
  not operational, to be thought again */
//   Sfield.time=Sfield.stream->time();
  Sfield.time=stream->time;
  Sfield.next=stream->next;
  
  const mesh_t
    *mesh=stream->mesh;
  if(stream->mesh==0)
    TRAP_ERR_RETURN(-1,1,"%s() called with Sfield.stream->mesh=%p\n",mesh);
  
  Sfield.units=strdup(stream->units);
  
/*-----------------------------------------------------------------------------
  then interpolate on structured one*/
  const double
    *data=stream->x;
  const list_t
    *list=stream->list;
  const grid_t
    *listGrid=stream->listGrid,
    *grid=Sfield.grid;
  const int
    discretisation=stream->descriptor->type;
  
  #pragma omp parallel
  {
  int i,j,m,e,e0=0;
  double lon,lat,*xm;
  
  #pragma omp for schedule(dynamic,1)
  for(j=0;j<grid->ny;j++){
    m=grid->nx*j;
    
    if(timeIsOld()){
      STDERR_BASE_LINE_FUNC("%s%3.0f%%%s",el,100.*j/grid->ny,cr);
      }
    
    for(i=0;i<grid->nx;i++,m++){
      xm=&Sfield.x[m];
      grid->xy(m,lon,lat);
      e=fe_detect(*mesh,mesh->triangles,*list,*listGrid,lon,lat,e0);
      fe_interpolate2D(*mesh,discretisation,data,stream->mask,lon,lat,e,xm);
      if(*xm==stream->mask)
        *xm=Sfield.mask;
      }
    
    }
  
  }
  STDERR_BASE_LINE_FUNC("%sdone.\n",el);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  SGfield_t<double> * InitDataStream(double time, const char *path, const char *convention, int format, grid_t *grid, double mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  mesh_t *mesh=new mesh_t;
  discretisation_t *descriptor;
  list_t *list=new list_t;
  grid_t *listGrid=new grid_t;
  int discretisation;
  int (* UGprocessor) (UGfield_t<double> & field, double t)=0;
  int (*identify)(const char *, const char *, int *, int);
  
  filestream_t<double> *zfilestream;
  sequence_t<UGfield_t<double> > *zUGsequence;
  UGfield_t<double>  *zUGstream;
  SGfield_t<double> *zSGstream;
  
  switch(format) {
    case ARCHIVE_FORMAT_BINARY:
      UGprocessor=stream_UGxBINARY;
      identify=identify_BINARY;
      break;
    case ARCHIVE_FORMAT_NETCDF:
      UGprocessor=stream_UGxNETCDF;
      identify=0;
      break;
    default:
      TRAP_ERR_EXIT(-1, "illegal format\n\n");
      break;
    }

  zfilestream=new filestream_t<double> (path,convention,"elevation",identify);
  status=zfilestream->finalize(time);
  
/*-----------------------------------------------------------------------------
  read mesh */
  switch(format) {
    case ARCHIVE_FORMAT_BINARY:{
      meta_archive_t info;
      status=clx_archivereadheader(zfilestream->filename.c_str(), &info);
      if(status != 0) TRAP_ERR_RETURN(0,1,"clx_archivereadheader() error %d\n",status);
      status=fe_list(&(info.mesh));
      if(status!=0) TRAP_ERR_RETURN(0,1,"fe_list() error %d\n",status);
      *mesh=info.mesh;
      }break;
    case ARCHIVE_FORMAT_NETCDF:
      status=fe_readmesh3d(zfilestream->filename.c_str(), mesh, 0);
      if(status!=0) TRAP_ERR_RETURN(0,1,"fe_readmesh3d() error %d\n",status);
      status=fe_geometry(mesh);
      if(status!=0) TRAP_ERR_RETURN(0,1,"fe_geometry() error %d\n",status);
      break;
    default:
      TRAP_ERR_EXIT(-1, "illegal format\n\n");
      break;
    }
  
/*-----------------------------------------------------------------------------
  initialise descriptor */
  discretisation=LGP1;
  
  status=discretisation_init(mesh,discretisation,0);
  descriptor=get_descriptor_address(*mesh,discretisation);
  
/*-----------------------------------------------------------------------------
  initialise UG field sequence */
  zUGsequence=new sequence_t<UGfield_t<double> >;
  status=zUGsequence->allocate(2);
  
  for(int k=0;k<zUGsequence->nframes;k++) {
/*-----------------------------------------------------------------------------
    initialise UG stream */
    zUGsequence->frames[k].descriptor=descriptor;
    zUGsequence->frames[k].processor=UGprocessor;
    zUGsequence->frames[k].allocate();
    zUGsequence->frames[k].stream=zfilestream;
    }
  status=zUGsequence->init(time);
  
/*-----------------------------------------------------------------------------
  initialise time interpolated UG stream */
  zUGstream=new UGfield_t<double>;
  
/*-----------------------------------------------------------------------------
  tricky bypass */
  zUGstream->stream=(datastream_t< double > *) zUGsequence;
//   zUGstream->stream=zUGsequence->buffer;
  
  fe_Allocate_and_CreateList(*mesh,listGrid,list);
  
  zUGstream->descriptor=descriptor;
  zUGstream->mesh=mesh;
  zUGstream->list=list;
  zUGstream->listGrid=listGrid;
  
  zUGstream->allocate();
  zUGstream->processor=stream_copy;
  status=zUGstream->check(time);
  
/*-----------------------------------------------------------------------------
  initialise space interpolated SG stream */
  zSGstream=new SGfield_t<double>;
  zSGstream->stream=zUGstream;
  zSGstream->grid=grid;
  zSGstream->mask=mask;

  zSGstream->allocate();
  zSGstream->processor=stream_SGxUG;
  
  return (zSGstream);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *ouput_name(date_t date, const char *path, const char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int hour,minute,mdum;
  float second;
  char *output=new char[1000];
  
  second=fmod(date.second,60.f);
  mdum=(date.second-second)/60.;
  minute=fmod(mdum,60.);
  hour=(mdum-minute)/60.;
  sprintf(output,"%s/%s-%04d-%02d-%02d_%02d-%02d-%04.1f.nc",path, rootname, date.year,date.month,date.day,hour,minute,second);
  
  return(output);
}


string cmd;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int interp_DEMs(const vector<string> & topoPaths,const grid_t & grid,double mask,double *depths,const string & output)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// interpolate DEMs
/**
\param[in,out] *depths array[grid.HSize()] of depths.
  All masked values will be interpolated.
  All unmasked values will be multiplied by -1.
*/
/*----------------------------------------------------------------------------*/
{
  const size_t tc=topoPaths.size();
  const int vc=2;
  
  STDERR_BASE_LINE("interpolating topography\n");
  int status,ti,vi,j;
  
  const string varnames[vc]={"z","bathymetry"};
  grid_t *topoGrids=0;
  
  if(tc>0)
    topoGrids=new grid_t[tc];
  
/*------------------------------------------------------------------------------
  load topographies */
  
  for(ti=0;ti<tc;ti++){
    const string *topoPathi=&topoPaths[ti];
    grid_t *topoi=&topoGrids[ti];
    
    poc_data_t<double> var;
    
    printf("topography from "+*topoPathi+" : trying");fflush(stdout);
    
    for(vi=0;vi<vc;vi++){
      printf(" "+varnames[vi]+" ...");fflush(stdout);
      status=var.init(*topoPathi,varnames[vi],-1);
      if(status==0)
        break;
      }
    NC_TRAP_ERROR(wexit,status,1,"poc_var_t::init() error");
    printf("ok.\n");
    
    status=poc_get_grid(*topoPathi,var.info,topoi,-1);
    NC_TRAP_ERROR(wexit,status,1,"poc_get_grid() error");
    
    status=map_completegridaxis(topoi,1);
    
    topoi->zmask=var.mask;
    
    swapval(var.data,topoi->z);
    var.destroy_data();
    }
  
/*------------------------------------------------------------------------------
  interpolate */
#define debug_topo_interp 0
#if debug_topo_interp == 0
  #pragma omp parallel
#endif
  {
  double lon,lat,zm;
  int i,ti;
  int64_t m,*accel=0;
  if(tc>0)
    accel=new int64_t[tc];
  int verbose;
  grid_t *topoi;
#if debug_topo_interp == 0
  verbose=0;
#endif
  
  #pragma omp for schedule(dynamic,1)
  for(j=0;j<grid.ny;j++){
    
    m=grid.nx*j;
    
    if(tc>0)
      aset(accel,tc,(int64_t)-1);
    
    if(timeIsOld())
      STDERR_BASE_LINE("%sinterpolating:%3.0f%%%s",el,100.*j/grid.ny,cr);
    
    for(i=0;i<grid.nx;i++,m++){
      
      if(depths[m]!=mask){
        depths[m]*=-1;
        continue;
        }
      
      grid.xy(i,j,lon,lat);
      
#if debug_topo_interp
      verbose=abs(lon-1.)<1e-3 and abs(lat-49.4)<1e-3;
      if(verbose)
        STDERR_BASE_LINE_FUNC("%s(%g;%g)\n",el,lon,lat);
#endif
      
      for(ti=0;ti<tc;ti++){
        topoi=&topoGrids[ti];
        zm=topoi->zmask;
        index_interpolation(*topoi,lon,lat,&accel[ti],topoi->z,topoi->zmask,&zm,0);
        if(zm!=topoi->zmask) break;
        }
      
      if(ti>=tc)
        zm=mask;
      
      depths[m]=zm;
      }
    }
  
  deletep(&accel);
  }
  
  if(tc<=0)
    return 0;
  
  for(ti=0;ti<tc;ti++){
    topoGrids[ti].free();
    }
  delete[]topoGrids;
  
  poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
  global<<poc_att_t("history",cmd);
  
  poc_var_t lav,lov,ev;
  global<<poc_att_t("description","DEM file");
  
  const poc_dim_t
    lad=poc_dim_t("latitude",grid.ny),
    lod=poc_dim_t("longitude",grid.nx);
  global<<lad<<lod;
  
  const nc_type
#define SENSIBLE_TYPES 0
#if SENSIBLE_TYPES == 0
    ele_nct=NC_DOUBLE;
#else
    ele_nct=NC_FLOAT;
#endif
  
  lav.init(lad,NC_DOUBLE);
  lav<<poc_att_t(_FillValue,mask);
  lav<<poc_att_t("units","degrees_north");
  
  lov.init(lod,NC_DOUBLE);
  lov<<poc_att_t(_FillValue,mask);
  lov<<poc_att_t("units","degrees_east");
  
  ev.init("elevation",ele_nct,"","",mask);
  ev.dimensions<<lad<<lod;
  
  global<<ev<<lav<<lov;
  
  printf("over-writing "+output+" ");fflush(stdout);
  unlink(output.c_str());
  
  status=poc_create(output,global,0,0);
  
  printf(":"+ev.name);fflush(stdout);
  status=poc_put_var(output,ev,depths);
  
  printf(","+lav.name);fflush(stdout);
  status=poc_put_var(output,lav,grid.y);
  
  printf(","+lov.name);fflush(stdout);
  status=poc_put_var(output,lov,grid.x);
  
  printf("\n");
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int extract(date_t date, const char *path, const char *convention, const char *outpath, int format, grid_t & grid,const vector<string> & topoPaths,const string & demPath,const vector<plg_t> & flood)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,m,n,status;
  SGfield_t<double> *zSGstream;
  char *output;
  double time;
  
  const double mask=32767.;
  double *topo=0;
  
  poc_global_t global("constructed around " __LINE_FILE_PACKAGE_REVISION);
  global<<poc_att_t("history",cmd);
  
  poc_var_t lav,lov,tv;
  poc_data_t<double> UGbathymetry,water_depth;
  poc_data_t<int8_t> landtype;
  
  for(i=0;i<2;i++){
    
    if(i>0)
      date.add(1./24.);
    
    time=cnes_time(date,'d');
    
    double time_data;
    
    if(i==0){
      zSGstream=InitDataStream(time, path, convention, format, &grid, mask);
      
/*------------------------------------------------------------------------------
      NetCDF header, from Denis.Blumstein@cnes.fr 26 Apr 2016 16:09 */
      global<<poc_att_t("description","Water depths");
      
      const poc_dim_t
        lad=poc_dim_t("latitude",grid.ny),
        lod=poc_dim_t("longitude",grid.nx),
        td=poc_dim_t("time",1);
      global<<lad<<lod<<td;
      
      /* TODO: test more sensible types */
#if SENSIBLE_TYPES == 0
      const nc_type
        wd_nct=NC_DOUBLE,
        lt_nct=NC_INT64,
        tv_nct=NC_INT64;/* especially this one, knowing the unit... */
      typedef int64_t tv_t;
#else
      const nc_type
        wd_nct=NC_FLOAT,
        lt_nct=NC_BYTE,
        tv_nct=NC_DOUBLE;
      typedef double tv_t;
#endif
      
      water_depth.info.init("water_depth",wd_nct,"","m",mask);
      water_depth.info.dimensions<<td<<lad<<lod;
      water_depth.init();
      
      lav.init(lad,NC_DOUBLE);
      lav<<poc_att_t(_FillValue,mask);
      lav<<poc_att_t("units","degrees_north");
      
      lov.init(lod,NC_DOUBLE);
      lov<<poc_att_t(_FillValue,mask);
      lov<<poc_att_t("units","degrees_east");
      
      tv.init(td,tv_nct);
      tv.name="simulation_time";
      tv<<poc_att_t(_FillValue,(tv_t)mask);/* TODO: check REAL necessity ... */
      tv<<poc_att_t("units","days");
      
      landtype.info.init("landtype",lt_nct,"","",mask);
      landtype.info.dimensions<<td<<lad<<lod;
      landtype.init();
      
      global<<water_depth.info<<lav<<lov<<tv<<landtype.info;
      
      n=grid.Hsize();
      if(grid.modeH!=1)
        TRAP_ERR_EXIT(ENOEXEC,"%s not coded for grid.modeH=%d!=1\n",__func__,grid.modeH);
      
/*------------------------------------------------------------------------------
      bathymetry */
      UGfield_t<double> *zUGstream=(UGfield_t<double> *) zSGstream->stream;
      sequence_t<UGfield_t<double> > *zUGsequence=(sequence_t<UGfield_t<double> > *)zUGstream->stream;
      filestream_t<double> *zfilestream=(filestream_t<double> *)zUGsequence->frames[0].stream;
      
      status=UGbathymetry.init(zfilestream->filename,"bathymetry");
      
      /* TODO: check discretisation is same for elevation and bathymetry */
      swapval(zUGstream->x,UGbathymetry.data);
      
      STDERR_BASE_LINE_FUNC("interpolating bathymetry\n");
      zSGstream->check(time);
      zSGstream->time=-INFINITY;
      
      swapval(zUGstream->x,UGbathymetry.data);
      topo=new double[n];
      swapval(zSGstream->x,topo);
      
      status=interp_DEMs(topoPaths,grid,mask,topo,demPath);
      }
    
/*------------------------------------------------------------------------------
    sea/river surface height */
    STDERR_BASE_LINE_FUNC("interpolating elevation\n");
    status=zSGstream->check(time);
    
    char *flooded=0;
    int count=0;
//     const int i0=min(73410800,n-1);
    const int i0=min(4588100,n-1);
    
    if(flood.size()>0){
      STDERR_BASE_LINE_FUNC("extrapolating elevation[%d]",i0);
      
      status=map_completegridaxis(&grid,2);
      fprintf(stderr,"(%g,%g): mask",grid.x[i0],grid.y[i0]);
      flooded=plg_TestInterior(grid.x,grid.y,n,flood,0,0);
      status=map_completegridaxis(&grid,1);
      
      fprintf(stderr,":%d, extrapolation",flooded[i0]);
      double
        *distance=new double[n],
        *extrapolated=new double[n];
      distance_to_nearest_unmasked(grid,zSGstream->x,mask,distance,extrapolated);
      
      fprintf(stderr,":%g->(%g)%g, compilation",zSGstream->x[i0],distance[i0],extrapolated[i0]);
      delete[]distance;
      
      for(m=0;m<n;m++){
        char *floodedm=&flooded[m];
        if(*floodedm==PLG_POINT_EXTERIOR) continue;
        double *x=&zSGstream->x[m];
        if(*x!=mask){
          *floodedm=PLG_POINT_EXTERIOR;
          continue;
          }
        *x=extrapolated[m];
        count++;
        }
      
      fprintf(stderr,":%d %g->%g, %d done\n",flooded[i0],zSGstream->x[i0],extrapolated[i0],count);
      
      delete[]extrapolated;
      }
    
    STDERR_BASE_LINE_FUNC("processing elevation\n");
    count=0;
    #pragma omp parallel for reduction(+:count)
    for(m=0;m<n;m++){
      const double *x=&zSGstream->x[m],*topom=&topo[m];
      bool masked=*x==mask or *topom==mask;
      double *depthm=&water_depth.data[m];
      
      if(masked){
        *depthm=mask;
        }
      else{
        *depthm= *x-topo[m];
        if(flooded!=0 and flooded[m]!=PLG_POINT_EXTERIOR){
          if(*depthm>0)
            count++;
          else{
            *depthm=mask;
            masked=true;
            }
          }
        }
      
      landtype.data[m]= not masked;
      }
    
    if(flooded!=0){
      STDERR_BASE_LINE_FUNC("%g+%g,%d extrapolated\n",topo[i0],water_depth.data[i0],count);
      delete[]flooded;
      }
    
/*------------------------------------------------------------------------------
    write data */
    
    output=ouput_name(date, outpath, "extracted");
    printf("over-writing %s ",output);fflush(stdout);
    unlink(output);
    
    status=poc_create(output,global,0,NC_NETCDF4);
    
    printf(":"+water_depth.info.name);fflush(stdout);
    status=water_depth.write_data(output);
    
    printf(","+lav.name);fflush(stdout);
    status=poc_put_var(output,lav,grid.y);
    printf(","+lov.name);fflush(stdout);
    status=poc_put_var(output,lov,grid.x);
    
    time_data=i/24.;
    printf(","+tv.name);fflush(stdout);
    status=poc_put_var(output,tv,&time_data);
    
    printf(","+landtype.info.name);fflush(stdout);
    status=landtype.write_data(output);
    
    printf("\n");
    delete[] output;
    }
  
  delete[] topo;
  water_depth.destroy_data();
  landtype.destroy_data();
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "  %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Convert TUGOm output to simuSWOT-compliant landscape files.\n"
    "  The output grid MUST be specified with EITHER -z OR -g option.\n"
    "\n"
    "OPTIONS\n"
    "  -h,--help  Show this help and exit.\n"
    "  -p  followed by input directory. Default: .\n"
    "  -c  followed by input convention. Default: analysis.YYYY-MM.nc\n"
    "  -b  followed by input topography. Several inputs are allowed, from highest to lowest priority.\n"
    "  -z  followed by output grid zone name\n"
    "  -g  followed by output grid file\n"
    "  -f  followed by flood limits polygone.\n"
    "  -d  followed by date as YYYY/MM/DD_HH:MN:SS[...]\n"
    "  -o  followed by output directory. Default: .\n"
    );
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  char *convention=NULL,*path=NULL,*keyword;
  string zone,floodPath;
  char *outpath=NULL;
  char *gridfile=NULL;
  vector<string> topoPaths;
  grid_t grid;
  vector<plg_t> flood;
  date_t date=NADate/*,start_date,end_date*/;
  int hours,minutes;
  float seconds;
  int format;
  
  cmd=fct_echo( argc, argv);

  n=1;

  while (n < argc) {
    keyword=argv[n];
    
    if( strcmp(keyword,"--help")==0 or strcmp(keyword,"-h")==0 ){
      print_help(argv[0]);
      exit(0);
      }
    
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
/*------------------------------------------------------------------------------
        input directory */
        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        input convention */
        case 'c' :
          convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        topography */
        case 'b' :
          topoPaths.push_back(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        region of extraction */
        case 'z' :
          zone= argv[n+1];
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        grid of extraction */
        case 'g' :
          gridfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        flood limits polygone */
        case 'f' :
          floodPath= argv[n+1];
          n++;
          n++;
          break;

/*------------------------------------------------------------------------------
        working date */
        case 'd' :
          keyword=argv[n+1];
          n++;
          n++;
          hours=minutes=seconds=0;
          sscanf(keyword,"%d%*[/-]%d%*[/-]%d%*c%2d:%2d:%f",&date.year,&date.month,&date.day,&hours,&minutes,&seconds);
          date.second=(hours*60.+minutes)*60.+seconds;
          break;

// /*------------------------------------------------------------------------
//         starting date */
//         case 's' :
//           s= strdup(argv[n+1]);
//           n++;
//           n++;
//           sscanf(s,"%d/%d",&start_date.month,&start_date.year);
//           free(s);
//           break;
// 
// /*------------------------------------------------------------------------
//         ending date */
//         case 'e' :
//           s= strdup(argv[n+1]);
//           n++;
//           n++;
//           sscanf(s,"%d/%d",&end_date.month,&end_date.year);
//           free(s);
//           break;

/*------------------------------------------------------------------------------
        directory for extraction files */
        case 'o' :
          outpath= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        STDOUT_BASE_LINE("unknown option %s\n",keyword);
        exit(-1);
        break;
      }
    }
  
/*------------------------------------------------------------------------------
  check arguments */
  status=0;
  
  if( isnad(date) ){
    printf("*** Please specify the date with -d ***\n");
    status=1;
    }
  
  if( (zone == "") == (gridfile == 0) ){
    printf("*** Please specify the grid with EITHER -z OR -g ***\n");
    status=1;
    }
  
  if(status!=0){
    print_help(argv[0]);
    exit(-1);
    }
  
/*------------------------------------------------------------------------------
  default arguments */
  if(zone != "") {
    bool identified;
    grid=get_zonegrid(zone.c_str(), & identified);
    if (not identified)
      TRAP_ERR_EXIT(-1,"get_zonegrid(\""+zone+"\",) error: wrong zone name\n");
    map_completegridaxis(&grid,1);
    }
  else
    zone=(string)strrchr0(gridfile,'/');
  
  if(gridfile != NULL) status=cdf_loadvargrid (gridfile,0,&grid);
  
  if(convention==NULL) {
    convention= strdup("analysis.YYYY-MM.nc");
    }
  
  if(path==NULL) {
    path= strdup(".");
    }
  
  if(outpath==NULL) {
    outpath = strdup(".");
    }
  
  if(floodPath!=""){
    status=plg_load(floodPath,PLG_FORMAT_UNKNOWN,flood);
    if(status!=0)
      TRAP_ERR_EXIT(status,"plg_load(\""+floodPath+"\",,) error (%d %s)\n",status,strerror(status));
    }
  
/*------------------------------------------------------------------------------
  compute */
  //fscanf(stdin,"%d",&status);
//   date_t date(2008,03,01);
  
  const string
    demPath=(string)outpath+"/grid-archives-"+zone+".nc";
  
  format=ARCHIVE_FORMAT_NETCDF;
  status=extract(date, path, convention, outpath, format, grid, topoPaths, demPath, flood);
  
  STDOUT_BASE_LINE("end of grid-archives ... \n");
  
  exit(0);
}
