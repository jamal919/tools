
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/

#include <stdio.h>
#include <string>
#include <unistd.h> /* for unlink */

#include "tools-structures.h"

#include "map.h"
#include "poc-time.h"
#include "datastream.h"
#include "tides.h" /* for d2s */
#include "poc-netcdf.hpp"
#include "poc-netcdf-data.hpp"
#include "imbrication.h"

#include "meteo.h"               /* for #ECMWF and #GRB_THROUGH_POC_NC */
#include "poc-grib.h"            /* for poc_grib_get_time() */


using namespace std;


/* Fonctions externes --------------------------------------------------------*/

extern int wrf_loadtime(const char *, int, double **, int *, date_t *, int);

/* Variables locales --------------------------------------------------------*/

#if TUGO
static int meteo_model=ECMWF;
#endif

/* Fonctions appellees depuis Fortran-----------------------------------*/

#ifdef _add_

extern int cdf_checkfile_();
#   pragma weak cdf_checkfile_ = cdf_checkfile

#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_stream_DecodeName_help(int cpu)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///prints convention help for tide_decode_atlasname()
/**
\param mode Default: 3. See tide_decode_atlasname()
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "CONVENTION\n"
    "  \"YYYY\" is replaced by the year.\n"
    "  \"MM\" is replaced by the month 01 to 12.\n"
    "  \"DD\" is replaced by the day of month 01 to 31.\n");
  if(cpu>=0){
    printf("  \"CCCC\" is replaced by the process number.\n");
    }
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_DecodeName(double t, int cpu,const string & filedir,const string & file_convention, char **filename, double *dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  char dummy[256], *pointer;
  date_t actual;
  double dt0=NAN;
  date_t cdf_reference;
  FILE *out;
  
#if TUGO
/*------------------------------------------------------------------------------
  filestream a weak point, time origine should be CNES or any other standard one */
  getmodeldate(t, &actual);
#else
/*------------------------------------------------------------------------------
  expect CNES time in seconds */
  actual=poctime_getdatecnes(t, 's');
#endif

/*------------------------------------------------------------------------------
  build the file name */

  (*filename) = new char[filedir.size()+1+file_convention.size()+1];
  
  if(filedir=="") {
    sprintf((*filename), "%s", file_convention.c_str());
    }
  else {
    sprintf((*filename), "%s/%s", filedir.c_str(), file_convention.c_str());
    }

  out = fopen(*filename, "r");
  status = (out == NULL);

  switch (status) {
    case 0:
/*------------------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      break;

    default:
/*------------------------------------------------------------------------------
      use format information*/
      
      pointer = strstr((*filename), "YYYY");
      if(pointer != NULL) {
        sprintf(dummy, "%04d", actual.year);
        strncpy(pointer, dummy, 4);
        dt0=d2s*366;
        }
      
      pointer = strstr((*filename), "MM");
      if(pointer != NULL) {
        sprintf(dummy, "%02d", actual.month);
        strncpy(pointer, dummy, 2);
        dt0=d2s*31;
        }
      pointer = strstr((*filename), "DD");
      if(pointer != NULL) {
        sprintf(dummy, "%02d", actual.day);
        strncpy(pointer, dummy, 2);
        dt0=d2s;
        }
      
      if(cpu<0)
        break;
      
      pointer = strstr((*filename), "CCCC");
      if(pointer != NULL) {
        sprintf(dummy, "%04d", cpu);
        strncpy(pointer, dummy, 4);
        }
    }
  
  if(dt!=0)
    *dt=dt0;
  
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_DecodeName(double t,const string & filedir,const string & file_convention, char **filename, double *dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=stream_DecodeName(t,-1,filedir,file_convention,filename,dt);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int stream_download_template(filestream_t<T> * filestream, double t, bool background,int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const char *tmpdir=getenv("TMPDIR");
  if(tmpdir==0) TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when TMPDIR not set\n",__func__);
  
  char *tmp;
  int status;
  size_t size;
  
  status = stream_DecodeName(t,tmpdir, filestream->convention, &tmp);
  
  status=get_file_size(tmp,&size);
  if(verbose>1)STDERR_BASE_LINE("cpu %d: get_file_size(\"%s\",) returned %d: %s\n",gCPU_ID,tmp,status,strerror(status));
  if(status==0) return EEXIST;
  
#if TUGO
  /* make sure you do not create the file before the others check for its existence */
  status=P_MPI_Barrier(MPI_COMM_WORLD);
  
  /* make sure only one CPU per node downloads the file */
  if(verbose>1)STDERR_BASE_LINE("cpu %d: creating %s\n",gCPU_ID,tmp);
  fclose(fopen(tmp,"w"));
  status=P_MPI_Barrier(MPI_COMM_WORLD);
  if(verbose>1)STDERR_BASE_LINE("cpu %d: unlink'ing %s\n",gCPU_ID,tmp);
  status=unlink(tmp);
  
  if(verbose>1)STDERR_BASE_LINE("cpu %d: unlink(\"%s\") returned %d errno=%d: %s\n",gCPU_ID,tmp,status,errno,strerror(errno));
#else
  status=!status;
#endif
  if(status==0){
    if(verbose>=0)STDERR_BASE_LINE("cpu %d: downloading "+filestream->filename+"\n",gCPU_ID);
    
    struct timeval before;
    gettimeofday(&before);
    
    char *src;
    status = stream_DecodeName(t,filestream->path, filestream->convention, &src);
    
    string cmd="scp ";
    if(verbose>1)
      cmd+="-v ";
    cmd+=(string)src+" "+tmpdir;
    
    if(background){
      /* waiting a little because filestream function will be called for the next variable,
      so soon enough for scp not to have enough time to write anything otherwise */
      cmd+=" & sleep 1";
      }
    
    if(verbose>0)STDERR_BASE_LINE("cpu %d:+ "+cmd+"\n",gCPU_ID);
    status=system(cmd.c_str());
    if(verbose>=0)STDERR_BASE_LINE("cpu %d:returned %d after %gs\n",gCPU_ID,status,difftime(before));
    if(status!=0 and verbose>=0)
      TRAP_ERR_EXIT(status,"download of %s to %s failed (%d %s)\n",src,tmpdir,status,strerror(status));
    }
  
#if TUGO
  P_MPI_Barrier(MPI_COMM_WORLD);
#endif
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int stream_cleanup_template(filestream_t<T> * filestream, double t,int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const char *tmpdir=getenv("TMPDIR");
  if(tmpdir==0) TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when TMPDIR not set\n",__func__);
  
  char *tmp;
  int status;
  
  status = stream_DecodeName(t,tmpdir, filestream->convention, &tmp);
  
  if(verbose>0)STDERR_BASE_LINE("cpu %d: unlink'ing %s\n",gCPU_ID,tmp);
  status=unlink(tmp);
  if(status!=0 and verbose>=0)STDERR_BASE_LINE("cpu %d: unlink(\"%s\") error (%d %s)\n",gCPU_ID,tmp,errno,strerror(errno));
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_stream_download_and_cleanup_help()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  printf("\n"
    "ENVIRONMENT\n"
    "  If the TMPDIR environment variable is set, "
    "copies of forcing files to TMPDIR will be progressively made and deleted as needed.\n"
    "  Clusters may increasingly (in 2018) set this environment variable themselves.\n"
    "  If the TMPDIR environment variable is not set:\n"
    " - you SHOULD set it to a fast (e.g., in 2018, solid state) drive "
    "if your forcing files are on a slow (e.g., still in 2018, overloaded NFS) drive\n"
    " - as this uses scp to copy, you MUST set it "
    "to an at least local (at best also fast, see above) drive "
    "if your forcing files are only accessible through e.g. an SSH server.\n"
    "\n");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> void stream_download_and_cleanup_template(filestream_t<T> * filestream, double t, double dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  download and cleanup forcing

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status __attribute__((unused)),verbose=1;
  
  status=stream_download_template(filestream,t,false,verbose);
  
  /* using 1.25 rather than 1.5 because 1.25*2 is not an integer
  and therefore will not miss a month */
  status=stream_download_template(filestream,t+1.25*dt,true,verbose);
  
  if(status==0){
    if(verbose>1)STDERR_BASE_LINE("cpu %d:cleaning up\n",gCPU_ID);
    stream_cleanup_template(filestream,t-1.25*dt,verbose);
    }
}
#if 0
extern void stream_download_and_cleanup(filestream_t<T> * filestream, double t, double dt);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void stream_download_and_cleanup(filestream_t<T> * filestream, double t, double dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  stream_download_and_cleanup_template(filestream,t,dt);
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void stream_download_and_cleanup(filestream_t<float> * filestream, double t, double dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  stream_download_and_cleanup_template(filestream,t,dt);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void stream_download_and_cleanup(filestream_t<double> * filestream, double t, double dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  stream_download_and_cleanup_template(filestream,t,dt);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void stream_download_and_cleanup(filestream_t<short> * filestream, double t, double dt)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  stream_download_and_cleanup_template(filestream,t,dt);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int stream_SGxNETCDF_template(SGfield_t<T> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  convenient processor to feed an structured field from archived netcdf/grib
  
  GRIB specificities:
    the GRIB files structure is more flexible but less prdictible than netcdf
    one. To avoid systematic parsing of files to collect data, the inquire
    function has been modified to collect meteo data information and time frame 
    vector at first GRIB file use

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  int frame, nframe, verbose=0;
  int varid, hour;
  date_t actual;
  double increment;
  double *cdftime=0;
  date_t cdf_reference;
  cdfvar_t cdfvar;
  filestream_t<T> *stream;
  bool debug=false;
  bool spokesman=(gCPU_ID==gCPU_MASTER);
  bool eof;
  
  double scale,offset;
  
  stream=(filestream_t<T> *) field.stream;

#if TUGO
/*------------------------------------------------------------------------------
  filestream a weak point, time origine should be CNES or any other standard one */
  getmodeldate(t, &actual);
#else
/*------------------------------------------------------------------------------
  expect CNES time in seconds */
  actual=poctime_getdatecnes(t, 's');
#endif

/*------------------------------------------------------------------------------
  build the meteo file name and open it*/
  status = stream->check(t);
  NC_TRAP_ERROR(return,status,1,"stream->check(%g) error",t);

/*------------------------------------------------------------------------------
  expect time to be returned in seconds*/
  status = poc_find_timevarid(stream->filename,&varid);
  
#if !TUGO
  const int meteo_model=ECMWF;
#endif
  
  switch (meteo_model) {
    case ECMWF:
    case WWIII:
      status= poc_gettime(stream->filename, &cdf_reference, &cdftime, &nframe, 0);
      break;
    case WRF:
      status =wrf_loadtime(stream->filename.c_str(), varid, &cdftime, &nframe, &cdf_reference,0);
      break;
    default:
      TRAP_ERR_RETURN(-1,1,"meteo_model type %d unknown\n",meteo_model);
      break;
    }
  if(status != 0)
    TRAP_ERR_RETURN(status,1,"[%d]_loadtime(\""+stream->filename+"\",%d,...) error %d\n",meteo_model,varid,status);
  
#if TUGO
  for(frame = 0; frame < nframe; frame++) {
    cdftime[frame] += elapsed(t_reference, cdf_reference);
    }
#endif
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  identify the meteo archive frame just before or at requested time t
  
  special case: convenient frame is the last one of the file, can be ok or
                some frames can be missing in the meteo file.
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  increment=get_timesampling(cdftime, -1.0, nframe, 0);
  
  frame = 0;
  eof=false;
  while(cdftime[frame] <= t) {
    if(cdftime[frame]  == t) {
      frame++;
      
      break;
      }
    frame++;
    if(frame == nframe) {
      if(cdftime[frame-1] +increment <t) {
        eof=true;
        printf("%s : t=%lf cdftime[frame-1]+increment=%lf\n", __func__, t,cdftime[frame-1]+increment);
        }
      break;
      }
    }

/*------------------------------------------------------------------------------
  frame is one step ahead the wanted time frame*/
  frame = frame - 1;
  
#if TUGO
  getmodeldate(cdftime[frame], &actual);
#else
  actual=poctime_getdatecnes(cdftime[frame], 's');
#endif
  hour = (int) round(actual.second / 3600.);
  
  if(frame==-1) {
    TRAP_ERR_EXIT(-1,"%s : file=%s, first time frame is posterior to requested time frame %2.2d/%2.2d/%4d %2.2dh\n",__func__,stream->filename.c_str(), actual.day,actual.month, actual.year, hour);
    }

  if(eof) {
    TRAP_ERR_EXIT(-1,"%s : file=%s, last time frame is anterior to requested time frame %2.2d/%2.2d/%4d %2.2dh\n",__func__,stream->filename.c_str(), actual.day,actual.month, actual.year, hour);
    }

  if(spokesman){
    printf("read meteo (netcdf/grib) at %2.2d/%2.2d/%4d %2.2dh from %s : %s\n", actual.day,actual.month, actual.year, hour, stream->filename.c_str(), stream->varname.c_str());
    if(debug) printf("frame number: %d, file time (from model start): %lfh, model time : %lfh \n", frame, cdftime[frame], t / 3600.0);
    fflush(stdout);
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  Read data at target time frame

  GRIB file may contain fields with evoluting resolution
  
  we provide a mecanism to deal with changes in field resolution
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  varid = stream->id;
  poc_var_t *var=&stream->glob.variables[varid];
  
  status=poc_decode_mask(*var,&scale,&offset,&field.mask);
  
#if TUGO  
  if(is_grib(stream->filename)) {
/*------------------------------------------------------------------------------
    check for possible change in GRIB variable size */
    size_t size;
    int length;
    status=poc_get_var_length(*var,&length,0,0);
    const poc_att_t *att=var->attributes.findP(POC_GRIB_FILE_SIZE_VATT_NAME);
    size=(*att)[frame];
    if(size!=length) {
      poc_dim_t dim=var->dimensions.back();
      const poc_att_t *att;
      att=var->attributes.findP(POC_GRIB_FILE_SIZE_VATT_NAME);
      dim.len=(*att)[frame];
      var->dimensions<<dim;
      field.grid->free();
      status=poc_get_grid(stream->filename,*var,field.grid,frame,verbose);
      field.allocate();
      }
    }
#endif
  
  status=poc_get_vara(stream->filename, *var, frame, field.x, verbose);
  
  if(status!=0){
    delete [] cdftime;
    TRAP_ERR_RETURN(status,verbose,"poc_get_vara(\""+stream->filename+"\",(\""+var->name+"\"),) error %d (see above)\n",status);
    }
  
  poc_scale_data(field.x, (int) field.grid->Hsize(), scale, offset, &field.mask, field.mask, verbose);
  
  if(spokesman and debug) printf("\n varname=%s, id=%d\n", var->name.c_str(), varid);
  if(status != 0) {
    TRAP_ERR_RETURN(-1,1,"error reading %s from "+stream->filename+"\n", var->name.c_str());
    }

  field.time=cdftime[frame];
  if(frame<nframe-1) {
    field.next=cdftime[frame+1];
    }
  else {
    field.next=cdftime[frame]+increment;
    }
  
  const poc_att_t *unitAtt=var->attributes.findP("units");
  field.units=poc_strdup(unitAtt->as_charp());
  
#if TUGO
  status=map_extendbuffer(*field.grid, &(field.x));
#endif
  
/*------------------------------------------------------------------------------
  apply land/sea mask*/
//   if(meteo_landsea!=NULL) {
//     for(k = 0; k < nx*ny; k++)
//       if(meteo_landsea[k] == '\1') {
//         field.x[k] = mask;
//         }
//     }

  delete [] cdftime;

  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_SGxNETCDF(SGfield_t<float> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=stream_SGxNETCDF_template(field, t);
  
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_SGxNETCDF(SGfield_t<short> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=stream_SGxNETCDF_template(field, t);
  
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_UGxSG(UGfield_t<float> & Ufield, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  interpolate a structured field on an unstructured mesh

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int n, status=0;
  double x,y;
  float  *data,z,mask;
  grid_t *grid;
  fieldstream_t<SGfield_t<float> > *stream;
  int masked=0;
  
  stream=(fieldstream_t< SGfield_t<float> > *) Ufield.stream;
  
/**-----------------------------------------------------------------------------
  update stream, i.e. structured field*/
  status=stream->data.acquire(t);
  if(status!=0) return(status);
  
  Ufield.time=Ufield.stream->time();
  Ufield.next=stream->data.next;
  
  mask=stream->data.mask;
  grid=stream->data.grid;
  data=stream->data.x;
  
  Ufield.units=poc_strdup(stream->data.units);
  
/**-----------------------------------------------------------------------------
  then interpolate on unstructured one*/
#if GRB_THROUGH_POC_NC
  if(grid->modeH==2)
    map_completegridaxis(grid,1);
  int64_t m=-1;
#endif
  for(n = 0; n < Ufield.descriptor->nnodes; n++) {
    x = Ufield.descriptor->nodes[n].lon;
    y = Ufield.descriptor->nodes[n].lat;
    
    x=map_recale(*grid,x);
#if GRB_THROUGH_POC_NC
    updatemin(&y,grid->ymax);
    updatemax(&y,grid->ymin);
    index_interpolation(*grid, x, y, &m, data, mask, &z, 0);
#else
    status = map_interpolation(*grid, data, mask, x, y, &z);
#endif
    if(status != 0 or z==mask) {
      z=1e+10;
      masked++;
      }
    Ufield.x[n] = z;
    }
  
/*------------------------------------------------------------------------------
  diffuse to avoid masked value */
#warning Done on 1st found SO NOT ALWAYS ON THE NEAREST UNMASKED: REVIEW THIS!!!
  while(masked!=0) {
    masked=0;
    for(n = 0; n < Ufield.descriptor->nnodes; n++) {
      if(Ufield.x[n]!=1e+10) continue;
      for(int k=0;k<Ufield.descriptor->nodes[n].nnghbs;k++) {
        int nn=Ufield.descriptor->nodes[n].nghbs[k];
         if(Ufield.x[nn]!=1e+10) {
          Ufield.x[n]=Ufield.x[nn];
          break;
          }
        }
      if(Ufield.x[n]==1e+10) masked++;
      }
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> int stream_UGxNETCDF_template(UGfield_t<T> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  read a UG field from a netcdf file; time frame is the one just preeceding time
  
  13/03/2016 : time handling needs to be strongly improved

------------------------------------------------------------------------------*/
{
  int status;
  int frame, nframes;
  filestream_t<T> *stream;
  poc_global_t global;
  poc_data_t<T> var;
  double *cdftime;
  date_t cdf_reference, origine;
  
  stream=(filestream_t<T> *) field.stream;
  
  if(stream==0) return(-1);
  
/*------------------------------------------------------------------------------
  build the input file name and open it*/
  status = stream->finalize(t,-1);
  if(status!=0) return status;
  
  status=poc_inq(stream->filename,&global,-1);
  
/*------------------------------------------------------------------------------
  read the full time vector */
  status= poc_gettime(stream->filename, &origine, &cdftime, &nframes, 0);
  if(status!=0) return status;
  
/**----------------------------------------------------------------------
  from TUGO, deprecated*/
//   for(frame = 0; frame < nframes; frame++) {
//     cdftime[frame] += elapsed(t_reference, cdf_reference);
//     cdftime[frame] /= 3600.;
//     }
// 
//   for(frame = 0; frame < nframes; frame++) {
//     cdftime[frame] /= 3600.*24.;
//     }

/**----------------------------------------------------------------------
  from TUGO, deprecated*/
//   frame = 0;
//   while(cdftime[frame] * 3600.0 <= t) {
//     frame++;
//     if(frame == nframes) break;
//     }

  frame = 0;
  while(cdftime[frame] <= t) {
    frame++;
    if(frame == nframes)
      break;
    }

/*------------------------------------------------------------------------------
  frame is one step ahead the wanted time frame*/
  frame = frame - 1;

  status=var.init(global,global.variables[stream->id].name);
  if(status!=0) return status;
  STDOUT_BASE_LINE("reading frame %d of "+stream->filename+" (%s)\n",frame, poctime_sdate_cnes(cdftime[frame],'s'));
  status=var.read_data(stream->filename, frame, 0, 1);
  if(status!=0) return status;
 
/**----------------------------------------------------------------------
  from TUGO, deprecated*/
//   field.time=cdftime[frame]* 3600.0;
//   if(frame<nframes-1) {
//     field.next=cdftime[frame+1]* 3600.0;
//     }
//   else {
//     field.next=(cdftime[frame]+cdftime[frame]-cdftime[frame-1])* 3600.0;
//     }
  
  field.time=cdftime[frame];
  if(frame<nframes-1) {
    field.next=cdftime[frame+1];
    }
  else {
    field.next=(cdftime[frame]+cdftime[frame]-cdftime[frame-1]);
    }
  
  valcpy(field.x,var.data,field.descriptor->nnodes);
  field.mask=var.mask;
  
  const poc_att_t *units=var.info.attributes.findP("units");
  if(units!=0)
    field.units=strdup(units->as_charp());
  
  delete[]cdftime;
  var.destroy_data();
  
  return (0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_UGxNETCDF(UGfield_t<float> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=stream_UGxNETCDF_template(field,t);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_UGxNETCDF(UGfield_t<double> & field, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=stream_UGxNETCDF_template(field,t);
  return status;
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> int filestream2UGfield_processor_template(UGfield_t<T> & UGfield,double modeltime)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,verbose=0;
  
/*-----------------------------------------------------------------------------
  update file stream */
  filestream_t<T> *filestream=(filestream_t<T> *)UGfield.stream;
  
  filestream->finalize(modeltime);
  
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  In the following sections of this function,
  structured grid and buffers :
    - are :
      + initialised just for this occurence
      + to be cleaned up at the end
    - even though they will be used at the next occurence
  However this kind of sequence (will) tend to be many and so
  saves the memory for the structured grid and buffers of the other parameters.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
  
/*-----------------------------------------------------------------------------
  init SG buffer (see note above) */
  poc_data_t<T> SGdata;
  
  status=SGdata.init(filestream->filename,filestream->varname);
  NC_TRAP_ERROR(wexit,status,verbose,"poc_data_t::init(\""+filestream->filename+"\",\""+filestream->varname+"\") error");
  
/*-----------------------------------------------------------------------------
  load SG grid (see note above) */
  grid_t grid;
  
  status=poc_get_grid(filestream->filename,SGdata.info,&grid);
  NC_TRAP_ERROR(wexit,status,verbose,"poc_get_grid(\""+filestream->filename+"\",(\""+filestream->varname+"\"),) error");
  set_grid_list(&grid,1,verbose);
  
/*-----------------------------------------------------------------------------
  get file times (see note above) */
  double *filetimes;
  size_t nframes;
  date_t origine;
  
  status=poc_gettime(filestream->filename,&origine,&filetimes,&nframes,0);
  NC_TRAP_ERROR(wexit,status,verbose,"poc_gettime(\""+filestream->filename+"\",...) error");
  
/*-----------------------------------------------------------------------------
  compute frame index and time increment (see note above) */
  int frame;
  double time;
  
#if TUGO
  time=cnestime_d(modeltime)*d2s;
#else
  time=modeltime;
#endif
  
  frame=vpos(time,filetimes,nframes);
  
  double dt=NAN;
  
  if(frame<0 or nframes<=frame)
    TRAP_ERR_EXIT(ENOEXEC,"programming error : frame=%d/%d\n",frame,nframes);
  else if(frame<nframes-1)
    dt=filetimes[frame+1]-filetimes[frame];
  else{
    /* NOTE: hoping time increment is constant */
    dt=filetimes[frame]-filetimes[frame-1];
    }
  
  delete[]filetimes;
  
/*-----------------------------------------------------------------------------
  load SG buffer */
  
  status=SGdata.read_data(filestream->filename,frame,0,1);
  NC_TRAP_ERROR(wexit,status,verbose,"poc_data_t::read_scaled_data(\""+filestream->filename+"\",%d) error",frame);
  
  T *SGbuffer=SGdata.data;
  const T mask=SGdata.mask;
  
/*-----------------------------------------------------------------------------
  interpolate */
  T *UGbuffer,z;
  int i;
  int64_t prior=-1;
  
  discretisation_t *descriptor=UGfield.descriptor;
  const int nnodes=descriptor->nnodes;
  node_t *node;
  double lat;
  
  UGbuffer=UGfield.x;
  UGfield.mask=mask;
  
  for(i=0;i<nnodes;i++){
    node=&descriptor->nodes[i];
    lat=node->lat;
    updatemax(&lat,grid.ymin);
    updatemin(&lat,grid.ymax);
    index_interpolation(grid,node->lon,lat,&prior,SGbuffer,mask,&z,verbose);
    UGbuffer[i]=z;
    }
  
/*-----------------------------------------------------------------------------
  record update */
  
  UGfield.time=modeltime;
  UGfield.next=UGfield.time+dt;
  
/*-----------------------------------------------------------------------------
  clean-up SG grid and data (see note above) */
  
  SGdata.destroy_data();
  grid.free();
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int filestream2UGfield_processor(UGfield_t<float> & UGfield,double time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=filestream2UGfield_processor_template(UGfield,time);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> sequence_t< UGfield_t<T> >  initialise_UG_sequence_template( double time,
                                                    filestream_t<T> *filestream,
                                                    int nframes,mesh_t *mesh, discretisation_t *descriptor,
                                                    int (* processor) (UGfield_t<T> & field, double t) )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  sequence_t< UGfield_t<T> > UGsequence;
  
/*-----------------------------------------------------------------------------
  initialise file stream */
  status=filestream->finalize(time);
  
/*-----------------------------------------------------------------------------
  initialise sequence stream */
  status=UGsequence.allocate(nframes);
  
  for(int k=0;k<UGsequence.nframes;k++) {
    UGfield_t<T> *field=new UGfield_t<T>(filestream,mesh,descriptor,processor);
    UGsequence.frames[k]=*field;
    }
  status=UGsequence.init(time);
  if(status != 0) {
    TRAP_ERR_EXIT(status,"sequence_t<>::init(%g) error\n",time);
    }
  
  return (UGsequence);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  sequence_t< UGfield_t<float> >  initialise_UG_sequence(double time,
                                              filestream_t<float> *filestream,
                                              int nframes,mesh_t *mesh, discretisation_t *descriptor,
                                              int (* processor) (UGfield_t<float> & field, double t) )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  sequence_t< UGfield_t<float> > sequence;
  sequence=initialise_UG_sequence_template(time,filestream,nframes,mesh,descriptor,processor);
  return sequence;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> sequence_t< UGfield_t<T> >  initialise_UG_sequence_template( double time,
                                                    string path, string convention, const string & name,
                                                    int nframes,mesh_t *mesh, discretisation_t *descriptor,
                                                    int (* processor) (UGfield_t<T> & field, double t),
                                                    int (*identify)(const char *, const char *, int *, int) )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  sequence_t< UGfield_t<T> > UGsequence;
  
  if(path==""){
    path=convention;
    const char *p,*c;
    p=path.c_str();
    c=strrchr0(p);
    convention=c;
    path.erase(c-p);
    }
  
/*-----------------------------------------------------------------------------
  initialise file stream */
  filestream_t<T> *filestream=new filestream_t<T> (path.c_str(),convention.c_str(),name.c_str(),identify);
  
/*-----------------------------------------------------------------------------
  initialise sequence stream */
  UGsequence=initialise_UG_sequence_template(time,filestream,nframes,mesh,descriptor,processor);
  
  return (UGsequence);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  sequence_t< UGfield_t<float> >  initialise_UG_sequence(double time,
                                              const string & path, const string & convention, const string & name,
                                              int nframes,mesh_t *mesh, discretisation_t *descriptor,
                                              int (* processor) (UGfield_t<float> & field, double t),
                                              int (*identify)(const char *, const char *, int *, int) )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  sequence_t< UGfield_t<float> > sequence;
  sequence=initialise_UG_sequence_template(time,path,convention,name,nframes,mesh,descriptor,processor,identify);
  return sequence;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_imbricationxUG(imbrication_t<double> & Ufield, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
//   double x,y;
  double  *data,mask __attribute__((unused));
  
  fieldstream_t< UGfield_t<double> > *stream;
  stream=(fieldstream_t< UGfield_t<double> > *) Ufield.stream;
  
/*------------------------------------------------------------------------------
  update stream, i.e. source unstructured field*/
  status=stream->data.acquire(t);
  if(status!=0) return(-1);
  
  Ufield.time=Ufield.stream->time();
  
  Ufield.next=stream->data.next;
  
  mask=stream->data.mask;
  data=stream->data.x;
  
  Ufield.units=poc_strdup(stream->data.units);
  
/*------------------------------------------------------------------------------
  then interpolate on model unstructured field*/
  for(n = 0; n < Ufield.nvalues; n++) {
//     x = Ufield.nesting->basics[n].lon;
//     y = Ufield.nesting->basics[n].lat;
    if(Ufield.nesting->basics[n].element==-1) {
      Ufield.x[n] = 0;
      }
    else {
      Ufield.x[n] = Ufield.nesting->basics[n].interpolate(data);
      }
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int stream_imbricationxUG_02(imbrication_t<double> & Ufield, double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n, status;
  double  *data,mask;
  
  UGfield_t<double> *stream;
  stream=( UGfield_t<double> *) Ufield.stream;
  
/*------------------------------------------------------------------------------
  update stream, i.e. source unstructured field*/
  status=stream->acquire(t);
  if(status!=0) return(-1);
  
  Ufield.time=stream->time;
  
  Ufield.next=stream->next;
  
  mask=stream->mask;
  data=stream->x;
  
  Ufield.units=poc_strdup(stream->units);
  
/**-----------------------------------------------------------------------------
  then interpolate on model unstructured field*/
  for(n = 0; n < Ufield.nvalues; n++) {
    Ufield.x[n] = Ufield.nesting->basics[n].interpolate(data, mask, 0.);
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 sequence_t< imbrication_t <double> >  initialise_imbrication_sequence(double time, string path, string convention, const char *name,
                                                                        int nnodes, int nframes, discretisation_t *descriptor,
                                                                        nesting_t<double> *nesting,
                                                                        int (* processor) (UGfield_t<double> & field, double t),
                                                                        int (*identify)(const char *, const char *, int *, int))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  wrapper

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
  sequence_t< imbrication_t <double> > OBCsequence;
  
/*-----------------------------------------------------------------------------
  initialise u-velocity file stream */
  filestream_t<double> *filestream=new filestream_t<double> (path.c_str(),convention.c_str(),name,identify);
  status=filestream->finalize(time);
  
/*-----------------------------------------------------------------------------
  initialise u-velocity UG stream */
  fieldstream_t< UGfield_t<double> > *fieldstream;
  fieldstream=new fieldstream_t< UGfield_t<double> >;
  
  fieldstream->data.init(filestream,0,descriptor,processor);
  
/*-----------------------------------------------------------------------------
  initialise u-velocity OBC stream */
  status=OBCsequence.allocate(nframes);
  
  for(int k=0;k<OBCsequence.nframes;k++) {
/**----------------------------------------------------------------------------
    OBCstream is a generic imbrication stream*/
    imbrication_t <double> *OBCstream= new imbrication_t <double> (nnodes, nesting, (datastream_t<double> *) fieldstream, stream_imbricationxUG);
    OBCsequence.frames[k]=*OBCstream;
    }
  status=OBCsequence.init(time);
  if(status != 0) {
    TRAP_ERR_EXIT(-1,"error \n");
    }
  return (OBCsequence);
}


