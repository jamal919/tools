
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief poc-netcdf prototypes

Declares functions defined in fe-netcdf.cpp netcdf-io.cpp netcdf-obsoletes.cpp netcdf-parse.cpp netcdf-standards.cpp netcdf-structured.cpp poc-netcdf01.cpp poc-netcdf02.cpp poc-netcdf03.cpp poc-netcdf04.cpp poc-netcdf05.cpp poc-netcdf06.cpp poc-netcdf07.cpp
*/
/*----------------------------------------------------------------------------*/

#if NETCDF_PROTO_H == 0
#define NETCDF_PROTO_H 1

#if TUGO
#include <config.h>
#endif

#include <netcdf.h>

#include "fe-classes.h" /* for mesh_t */
#include "poc-time-classes.h" /* for date_t */

#include "netcdf-classes.h"
#include "poc-netcdf-assertions.h" /* for nc_check_error */

#if TUGO
#else
#include "tools-structures.h" /* for spectrum_t */
#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  list of time variable names

  see poc_gettime(int,date_t*,double**,size_t*,cdfvar_t*) and isT() for its uses.

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

static const char *timeNames[] __attribute__((unused)) ={
  "time",                         //for pretty much everything
  "time_counter","time_instant",  //for MERCATOR
  "simulation_time",              //for SWOT
  "scrum_time"                    //for ROMS
  ,"t","nT",
  "month",
  "frequency","W","ntide"
  };

#define LASTFRAME            -1
#define NOFRAME              -2

#define LASTLEVEL            -1
#define NOLEVEL              -2

#include "map.h" /* for grid_t */

//from netcdf-vio.cpp

extern int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],char *z);
extern int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],signed char *z);
extern int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],short *z);
extern int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],int *z);
extern int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],float *z);
extern int poc_get_vara(int ncid, int id, const size_t start[], const size_t count[],double *z);

extern int poc_get_var(int ncid, int id, char   *z);
extern int poc_get_var(int ncid, int id, int8_t *z);
extern int poc_get_var(int ncid, int id, short  *z);
extern int poc_get_var(int ncid, int id, float  *z);
extern int poc_get_var(int ncid, int id, double *z);

extern int poc_get_var(int ncid,const string & var, double *z, int verbose=0);
extern int poc_get_var(int ncid,const string & var, float *z, int verbose=0);
extern int poc_get_var(int ncid,const string & var, short *z, int verbose=0);

extern int poc_get_var(const string & filename,const string & var, double *z, int verbose=0);
extern int poc_get_var(const string & filename,const string & var, float *z, int verbose=0);
extern int poc_get_var(const string & filename,const string & var, short *z, int verbose=0);


extern int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count, const char   *buffer);
extern int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count, const signed char *buffer);
extern int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count, const short  *buffer);
extern int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count, const int    *buffer);
extern int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count, const float  *buffer);
extern int poc_put_vara(int ncid, int varid, const size_t *start, const size_t *count, const double *buffer);

extern int poc_put_var(int ncid, int varid, const char   *buffer);
extern int poc_put_var(int ncid, int varid, const signed char *buffer);
extern int poc_put_var(int ncid, int varid, const short  *buffer);
extern int poc_put_var(int ncid, int varid, const int    *buffer);
extern int poc_put_var(int ncid, int varid, const float  *buffer);
extern int poc_put_var(int ncid, int varid, const double *buffer);

extern int poc_put_var(const char *filename, int varid, const short *out);
extern int poc_put_var(const char *filename, int varid, const int *out);
extern int poc_put_var(const char *filename, int varid, const float *out);
extern int poc_put_var(const char *filename, int varid, const double *out);


/*******************************************************************************
       MORE OR LESS OLD STUFF
*******************************************************************************/


//from netcdf-utilities.cpp
extern int cdf_info(const char* filename);

extern int cdf_attinfo(const int ncid, const int var, const int att, cdfatt_t *info, int verbose=0);
extern void cdf_print_attinfo(const cdfatt_t info);

extern int cdf_varinfo(const int ncid, const int var, cdfvar_t *info, int verbose=0);
extern int cdf_varinfo(const int ncid, const char *varname, cdfvar_t *info, int verbose);

extern int cdf_varinfo(const char *filename, const char *varname, cdfvar_t *info, int verbose=0);
extern int cdf_varinfo(const char *filename, int var, cdfvar_t *info);
extern int cdf_varinfo(const cdfgbl_t &global, const char *varname, cdfvar_t *info);
extern void cdf_print_varinfo(const cdfvar_t &info);

extern int cdf_diminfo(const int ncid, const int dimid, cdfdim_t *info, int verbose=0);
extern void cdf_print_diminfo(const cdfdim_t info);

extern int cdf_globalinfo(int ncid,cdfgbl_t *info,int verbose);
extern int cdf_globalinfo(const char *filename,cdfgbl_t *info,int verbose=0);

//#include <vector>
//extern int parse_filelist(vector< string > filelist, char** varnames, int nvars, date_t start, date_t final, double* tmin=NULL, double* tmax=NULL);
extern int cdf_identify(const char *filename,const char *varname);
extern int cdf_identify(cdfgbl_t g_info, const char *name);
extern int cdf_varid(const char *filename, const char *varname);
extern int cdf_varid(int ncid, const char *varname);
extern int cdf_identify_dimension(cdfgbl_t g_info, const char *name);


//from poc-netcdf01.cpp

extern int cdf_createvariable(const char *filename,cdfvar_t *variable);
//gloabl etribute
extern int cdf_put_global_att_float(char *filename,char *name,nc_type xtype,size_t len, float *fp);

//from poc-netcdf02.c
extern int cdf_loadvargrid_3d (const char* filename, int vdata, grid_t *grid);
extern int cdf_loadvargrid_2d (const char* filename, int vdata, grid_t *grid,int verbose=0);

extern int poc_getgrid3d  (const char* filename, cdfgbl_t global, cdfvar_t info, grid_t *grid, int *);

extern int poc_getgrid2d  (const char* filename, cdfgbl_t global, cdfvar_t info, grid_t *grid);
extern int poc_getgrid2d(const char *input,const char *variable, grid_t *grid,int verbose=0);

extern int poc_UGdecode_discretisation(cdfvar_t variable, cdfgbl_t global, UGdecoded_t *decoded);
extern int poc_UGdecode_axis(cdfvar_t variable, cdfgbl_t global, UGdecoded_t *decoded);

//from poc-netcdf03.c
extern int init_ncvariable(cdfvar_t *variable);
//extern int free_ncvariable(cdfvar_t *variable) __attribute__((safe_deprecated("use cdfvar_t::destroy() instead")));
//extern int delete_ncvariable(cdfvar_t *variable);
extern int create_ncfile3d(char *filename, grid_t grid);
extern int poc_createfile(const char *filename, int mode=-1);
extern int poc_createfile(const char *filename, const cdfgbl_t &global);

extern int poc_standardaxis(const char* vname, const char* gname, const char* units, const char* sname, const char* lname, double vmin, double vmax, char** dname, cdfvar_t* variable);
extern int poc_standardaxis_xy(const char* vname, const char* gname, const char* units, const char* sname, const char* lname, double vmin, double vmax, char** dname, cdfvar_t* variable);
extern int poc_standardaxis_xieta(char *vname,char *gname, char *units,char *sname,char *lname,double vmin,double vmax, char **dname, cdfvar_t *variable);
extern int poc_standardaxis_z(const char *name,float mask,const char *units,const char *standardname,const char *longname, char **dname, cdfvar_t *variable,const char *,const char *);
extern int poc_standardaxis_s(char *name,float mask,char *units,char *standardname,char *longname, char **dname, cdfvar_t *variable, char *gridname);
extern int poc_standardaxis_incidence(const char *name, const char *standardname,const char *longname, const char **dname, cdfvar_t *variable);

extern char *poc_getdate(date_t actual);
extern cdfvar_t poc_standardtime(const char *vname,const char *dname,const char *units,date_t origin);
extern void poc_standardvariable_xy(cdfvar_t *variable,const char *name, float mask,const char *units,float scale,float offset, const char *standardname,const char *longname,const char *shortname,pocgrd_t grid);

extern int poc_standardvariable_xyt(const char *name, float mask, const char *units, float scale,  float offset,  const char *standardname,const char *longname,const char *shortname, pocgrd_t grid, cdfvar_t &variable);
extern int poc_standardvariable_xyt(const char *name, short mask, const char *units, float scale,  float offset,  const char *standardname,const char *longname,const char *shortname, pocgrd_t grid, cdfvar_t &variable);
extern int poc_standardvariable_xyt(const char *name, double mask,const char *units, double scale, double offset, const char *standardname,const char *longname,const char *shortname, pocgrd_t grid, cdfvar_t &variable);

extern cdfvar_t poc_standardvariable_xyw(char *name, float mask,char *units,float scale,float offset,char *standardname,char *longname,char *shortname, pocgrd_t grid);
extern cdfvar_t poc_standardvariable_xywt(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname,pocgrd_t grid);
extern int poc_standardvariable_xyz(const char* name, float mask, const char* units, float scale, float offset, const char* standardname, const char* longname, const char* shortname, pocgrd_t grid, cdfvar_t &);
extern cdfvar_t poc_standardvariable_xyzw(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname,pocgrd_t grid);
extern int poc_standardvariable_xyzt(const char *name, float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname, pocgrd_t grid, cdfvar_t &);
extern cdfvar_t poc_standardvariable_xyzwt(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname,pocgrd_t grid);
extern cdfvar_t poc_standardvariable_xieta(char *name, float mask,char *units,float scale,float offset,char *standardname,char *longname,char *shortname, pocgrd_t grid);
extern cdfvar_t poc_standardvariable_xietat(char *name, float mask,char *units,float scale,float offset,char *standardname,char *longname,char *shortname, pocgrd_t grid);
extern cdfvar_t poc_standardvariable_xietaw(char *name, float mask,char *units,float scale,float offset,char *standardname,char *longname,char *shortname, pocgrd_t grid);
extern cdfvar_t poc_standardvariable_xietawt(char *name, float mask,char *units,float scale,float offset, char *standardname,char *longname,char *shortname, pocgrd_t grid);
extern cdfvar_t poc_standardvariable_xietas(char *name, float mask,char *units,float scale,float offset,char *standardname,char *longname,char *shortnamee, pocgrd_t grid);
extern cdfvar_t poc_standardvariable_xietasw(char *name, float mask,char *units,float scale,float offset,char *standardname,char *longname,char *shortnamee, pocgrd_t grid);
extern cdfvar_t poc_standardvariable_xietast(char *name, float mask,char *units,float scale,float offset,
                                          char *standardname,char *longname,char *shortname, pocgrd_t grid);
extern cdfvar_t poc_standardvariable_xietaswt(char *name, float mask,char *units,float scale,float offset,
                                          char *standardname,char *longname,char *shortname, pocgrd_t grid);
extern cdfvar_t poc_standardvariable_ift(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname);
extern int create_ncvariable(const char *filename,cdfvar_t *variable);
extern int create_ncvariable(int ncid,cdfvar_t *variable, int verbose=0);

extern void add_att(cdfvar_t *v,cdfatt_t a);

extern int create_ncvariable_indfile(char *filename,cdfvar_t *variable,char *filegrid);

extern int poc_sphericalgrid_xy   (const char *input, const char *name,const grid_t & grid, pocgrd_t *cdfgrid);
extern int poc_sphericalgrid_xy   (const char *input, pocgrd_t *cdfgrid);

extern int poc_sphericalgrid_xyz  (const char *input, const char *hlabel, const char *vlabel, grid_t grid, pocgrd_t *cdfgrid);
extern int poc_inq_grid(const char *input, const char *hlabel, const char *vlabel, pocgrd_t *cdfgrid);

extern int poc_sphericalgrid_xyt (const char *input, const char *name, grid_t grid, pocgrd_t *cdfgrid);

extern int poc_sphericalgrid_xyzt(const char *input, const char *hlabel, const char *vlabel, grid_t grid, pocgrd_t *cdfgrid);
extern int poc_sphericalgrid_xyzwt(const char *input, const char *name, int nw, grid_t grid, pocgrd_t *cdfgrid);

extern int poc_sphericalgrid_xietast(char *input, char *name, grid_t grid, pocgrd_t *cdfgrid);
extern int poc_sphericalgrid_xietaswt(char *input, char *name, int nw, grid_t grid, pocgrd_t *cdfgrid);

//from poc-netcdf04.c
extern int cdf_loadvargrid (const char* filename,int vdata, grid_t *grid);
extern int cdf_getvariable (int ncid,int *nvar, char ***name);

extern int cdf_loadvar_r1_3d (const char* filename, int vdata, int t, grid_t grid, float *buf, float *mask,variable_t *info);
extern int cdf_loadvar_r1_2d (const char* filename, int v, int k, int t, grid_t grid, int n, float *buf, float *mask,variable_t *info);

extern void cdff_saver1 (char* filename, int *k,int *d,int *t, int *tmax,grid_t *grid, int *n, float *buffer, float *mask, int *status);

//from poc-netcdf05.c
extern int create_header(char *output, grid_t grid, spectrum_t spectrum);
extern int create_headershort(char *output, grid_t grid, spectrum_t spectrum);
extern int nc_writeaxis(int ncid, grid_t grid) ;
extern int nc_writeaxis_3d(int ncid, grid_t grid) ;
extern int nc_writespectrum(int ncid, spectrum_t spectrum) ;
extern int nc_write_depth(int ncid, grid_t grid, int var, float *z, float mask);
extern int nc_write_r1(int ncid, grid_t grid, int var, int wave, float *z, float mask) ;
extern int nc_write_c1(int ncid, grid_t grid, int var, int wave, fcomplex *z) ;

#define ELEVATION 0
#define CURRENTS  1
#define PRESSURE  2
#define IBD       3
#define create_ncfile2d_OPTIONS 4

extern int create_ncfile2d(char *filename,size_t x_len, size_t y_len, grid_t grid,int options[create_ncfile2d_OPTIONS]);

extern bool isT(const char *name);
extern bool isT(const cdfdim_t &dim);

extern int poc_write(const char *filename, cdfvar_t var, size_t *start,size_t *count,char   *buffer);
extern int poc_write(const char *filename, cdfvar_t var, size_t *start,size_t *count,short  *buffer);
extern int poc_write(const char *filename, cdfvar_t var, size_t *start,size_t *count,float  *buffer);
extern int poc_write(const char *filename, cdfvar_t var, size_t *start,size_t *count,double *buffer);


extern int write_ncfile2dframe(char *filename, grid_t grid, int frame, int id, float *buffer);

extern int poc_write_xy(const char *filename, grid_t grid, int id,  short *buffer);
extern int poc_write_xy(const char *filename, grid_t grid, int id,    int *buffer);
extern int poc_write_xy(const char *filename, grid_t grid, int id,  float *buffer);
extern int poc_write_xy(const char *filename, grid_t grid, int id, double *buffer);

extern int poc_write_xyz(const char *filename, grid_t grid, int id, float *buffer);
extern int poc_write_xyz(const char *filename, grid_t grid, int id, double *buffer);

extern int poc_write_xyt(const char *filename, grid_t grid, int frame, int id, float *buffer);
extern int poc_write_xyt(const char *filename, grid_t grid, int frame, int id, double *buffer);

extern int poc_write_xyzt(const char *filename, grid_t grid, int frame, int id, float *buffer);
extern int poc_write_xyzt(const char *filename, grid_t grid, int frame, int id, double *buffer);

extern int write_ncfile2d(char *filename, grid_t grid, int ndim, int *tstart, int id, float *buffer);
extern int write_ncfile2d_obsolete(char *filename, grid_t grid, int frame, int id, float *buffer);
extern int write_ncfile3d(char *filename, grid_t grid, int ndim, int *tstart, int id, float *buffer);
extern int write_ncfile3d_obsolete(char *filename, grid_t grid, int frame, int id, float *buffer);

extern int poc_writetime(const char *filename, int frame, int id, double time);
extern int poc_writetime(const char *filename, int frame, const char *name, double time);
extern int poc_settime(const char *filename, int id, double *time);

extern int save_SG(const char *output,grid_t grid,float *buffer,float mask,const char *varname,const char *unit,const char *standard,int create);
extern int save_SG(const char *output,grid_t grid,fcomplex *tide,fcomplex cmask,const char **varname,const char **unit,const char **standard,int create);

extern int save_SGXY(const char *output, const grid_t & grid,    int *buffer,    int mask, const char *varname,const char *unit, const char *standard, int create, int write_grid=1);
extern int save_SGXY(const char *output, const grid_t & grid,  float *buffer,  float mask, const char *varname,const char *unit, const char *standard, int create, int write_grid=1);
extern int save_SGXY(const char *output, const grid_t & grid, double *buffer, double mask, const char *varname,const char *unit, const char *standard, int create, int write_grid=1);

extern int read_SGatlas(const char *filename,const char **varnames,grid_t & grid,complex<float> * tide,complex<float> *mask);
extern int read_SGatlas(const char *filename,const char **varnames,grid_t & grid,complex<float> * tide,complex<float> *mask);

//from netcdf-time.cpp
extern int poc_gettime(int ncid, int tvid, date_t *origine, double **time, size_t *nframes);
extern int poc_gettime(int ncid, date_t *origine, double **time, size_t *nframes, cdfvar_t *timevar=NULL);
extern int poc_gettime(const char *filename, date_t *origine, double **time, size_t *nframes, cdfvar_t *timevar=NULL);
extern int poc_gettime(int file, date_t *origine, double **time, int *nframes, cdfvar_t *timevar=NULL);
extern int poc_gettime(const char *file, date_t *origine, double **time, int *nframes, cdfvar_t *timevar);
extern int poc_gettime(int ncid, cdfvar_t variable, cdfgbl_t global, date_t *origine, double **time, size_t *nframes);
extern int poc_gettime(const char *filename, cdfvar_t variable, cdfgbl_t global, date_t *origine, double **time, size_t *nframes);

extern int poc_timefilterlist(vector<string> *pathlist,double *start,double *final,double **ts=NULL,double origd=NAN);
extern int getSkippedFrames(size_t nt,double *ts,int handleRepeatedFrames,bool **skipFrame);

extern int poc_getvar2d(const char* filename, int varid, int frame, signed char *buf, signed char *mask, cdfvar_t info);
extern int poc_getvar2d(const char* filename, int varid, int frame, short *buf, short *mask, cdfvar_t info);
extern int poc_getvar2d(const char* filename, int varid, int frame, float *buf, float *mask, cdfvar_t info);
extern int poc_getvar2d(const char* filename, int varid, int frame, double *buf, double *mask, cdfvar_t info);


extern int poc_getvar3d(const char* filename, int varid, int frame, float *buf, float *mask,cdfvar_t info);

//from poc-netcdf07.c
extern int poc_getshort_nt(const char *filename, mesh_t mesh, int frame, int id, float *z, float *mask);

//extern int poc_put_UG2D(const char *filename, mesh_t & mesh, int id, float *z);
extern int poc_put_UG2D(const char *filename, const mesh_t &, int, const int *);
extern int poc_put_UG2D(const char *filename, const mesh_t &, int, const float  *);
extern int poc_put_UG2D(const char *filename, const mesh_t &, int, const double *);


//extern int poc_put_UG3D(const char *filename, mesh_t, int, int, double *);

extern int poc_put_UG3D(const char *filename, const mesh_t &, int, int, const int *);
extern int poc_put_UG3D(const char *filename, const mesh_t &, int, int, const float  *);
extern int poc_put_UG3D(const char *filename, const mesh_t &, int, int, const double *);

//extern int poc_put_UG3D(const char *filename, mesh_t &, int, int, double *);
extern int poc_put_UG3D(const char *filename, mesh_t &, int, double **);

extern int poc_get_att(int ncid, int id , const char* name, short *z);
extern int poc_get_att(int ncid, int id , const char* name, int *z);
extern int poc_get_att(int ncid, int id , const char* name, float *z);
extern int poc_get_att(int ncid, int id , const char* name, double *z);

extern int poc_seqget_UG3D(const char *filename, mesh_t &, int, int, double *);

extern int poc_get_UG3D(int ncid, int frame, int id, int    *z);
extern int poc_get_UG3D(int ncid, int frame, int id, float  *z);
extern int poc_get_UG3D(int ncid, int frame, int id, double *z);

extern int poc_get_UG3D(int ncid, int frame, const char *varname, int    *z, int verbose=0);
extern int poc_get_UG3D(int ncid, int frame, const char *varname, double *z, int verbose=0);

extern int poc_get_UG3D(int ncid, int frame, const char *varname, complex<double> *z, double *aBuf,double *GBuf);

extern int poc_get_UG3D(const char* filename, int frame, const char* varname, int    *z, int verbose=0);
extern int poc_get_UG3D(const char *filename, int frame, const char *varname, float  *z, int verbose=0);
extern int poc_get_UG3D(const char *filename, int frame, const char *varname, double *z, int verbose=0);

extern int poc_get_UG3D(const char *filename, int frame, const char *amplitude, const char *phase, complex<float> *z);
extern int poc_get_UG3D(const char *filename, int frame, const char *amplitude, const char *phase, complex<double> *z);

extern int poc_get_UG3D(const char *filename, int frame, int level, const char *amplitude, const char *phase, complex<float> *z);
extern int poc_get_UG3D(const char *filename, int frame, int level, const char *amplitude, const char *phase, complex<double> *z);

extern int poc_put_UG4D(const char *filename, mesh_t, int, int, double **);
extern int poc_put_UG4D(const char *filename, mesh_t, int, int, float **);
extern int poc_put_UG4D(const char *filename, mesh_t, int, int, float *);

extern int poc_get_UG4D(const char *filename, mesh_t &, int, int, double **);
extern int poc_get_UG4D(int ncid, mesh_t &, int, int, double **);

extern  int poc_get_UG4D(const char *filename, int frame, const char *A, const char *G, complex<double> **z);

extern int poc_getLast(int ncid, int frame, int frequency, int id, short  *z);
extern int poc_getLast(int ncid, int frame, int frequency, int id, int    *z);
extern int poc_getLast(int ncid, int frame, int frequency, int id, float  *z);
extern int poc_getLast(int ncid, int frame, int frequency, int id, double *z);

extern int poc_putLast(int ncid, int frame, int frequency, int id, short  *z);
extern int poc_putLast(int ncid, int frame, int frequency, int id, int    *z);
extern int poc_putLast(int ncid, int frame, int frequency, int id, float  *z);
extern int poc_putLast(int ncid, int frame, int frequency, int id, double *z);

extern int poc_putLast(const char *filename, int frame, int frequency, int id, short  *z);
extern int poc_putLast(const char *filename, int frame, int frequency, int id, int    *z);
extern int poc_putLast(const char *filename, int frame, int frequency, int id, float  *z);
extern int poc_putLast(const char *filename, int frame, int frequency, int id, double *z);


extern cdfvar_t poc_shortvariable_nt(char *, float, char *, float, float,char *);
extern cdfvar_t poc_floatvariable_nt(char *, float, char *, float, float,char *);
extern cdfvar_t poc_floatvariable_nz(const char *name, float mask,const char *units,float scale,float offset,const char *standardname);

extern cdfvar_t poc_variable_UG2D(const char *, float,  const char *,float,  float, const char *,const char *);
extern cdfvar_t poc_variable_UG2D(const char *, double, const char *,double, double,const char *,const char *);
///extern cdfvar_t poc_variable_UG3D(const char *, double, const char *,double, double,const char *,const char *,const char *);

extern cdfvar_t poc_variable_UG3D(const char *, int,    const char *, int,    int,   const char *,const char *,const char *);
extern cdfvar_t poc_variable_UG3D(const char *, float,  const char *, float,  float, const char *,const char *,const char *);
extern cdfvar_t poc_variable_UG3D(const char *, double, const char *, double, double,const char *,const char *,const char *);

extern cdfvar_t poc_variable_UG3D(const char *, float,  const char *, float,  float, const char *,const char *,int);
extern cdfvar_t poc_variable_UG3D(const char *, double, const char *, double, double,const char *,const char *,int);

extern int      poc_variable_UG3D(const char *, int, const char *, int,    const char *, int,    int,   const char *, const char *, int);
extern int      poc_variable_UG3D(const char *, int, const char *, double, const char *, double, double,const char *, const char *, int);

extern cdfvar_t poc_variable_UG4D(const char *, float,  const char *, float,  float, const char *, const char *, const char *,const char *);
extern cdfvar_t poc_variable_UG4D(const char *, double, const char *,double, double, const char *, const char *, const char *,const char *);
extern cdfvar_t poc_variable_UG4D(const char *, double, const char *, double, double,const char *, const char *, int,         int);

extern int poc_get_timeorigin (const char* filename, int vtime, date_t *reference);

//extern int fe_ncinfo (char* filename, double *time);
extern int fe_extract(int ncid, int varid, size_t *start,size_t *count, float *buffer) ;
extern int fe_read2d(int ncid,int varid, int layer, int frame, float **buffer, float *mask) ;

extern int fe_read3d(int ncid,int varid, int frame, float **buffer, float *mask) ;
extern int fe_read3d(const char *filename, int varid, int frame, float **data, float *mask) ;


extern int poc_getbuffer(int ncid,int varid, float *buffer);
extern int poc_getbuffer(int ncid,const char*varname, float *buffer);

extern int poc_decode_axis(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded);
extern int decoded_nvalues(decoded_t decoded, int *ndims=NULL);
extern int poc_decode_associates(cdfvar_t variable, cdfgbl_t global, decoded_t *decoded, int initialize, bool debug) __attribute__((safe_deprecated("not sure whether this works with ROMS")));

extern int poc_decode_mask(const cdfvar_t &variable, signed char *scale, signed char *offset, signed char *spec);
extern int poc_decode_mask(const cdfvar_t &variable, char *scale, char *offset, char *spec);
extern int poc_decode_mask(const cdfvar_t &variable, short int *scale, short int *offset, short int *spec);
extern int poc_decode_mask(const cdfvar_t &variable, int *scale, int *offset, int *spec);
extern int poc_decode_mask(const cdfvar_t &variable, float *scale, float *offset, float *spec);
extern int poc_decode_mask(const cdfvar_t &variable, double *scale, double *offset, double *spec);
extern int poc_decode_mask(const cdfvar_t &variable, decoded_t *decoded);

extern int poc_decode_units(cdfvar_t variable, decoded_t *decoded, double *factor);
extern int poc_decode_names(cdfvar_t variable, decoded_t *decoded);
extern int poc_decode_names(cdfgbl_t global, decoded_t *decoded);


/*******************************************************************************
       This template is already in tugo
*******************************************************************************/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_getvar2d(const char* filename, T *buf, T *mask, cdfvar_t info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_getvar2d()
/*----------------------------------------------------------------------------*/
{
  int status;
  
  ///get the 1st frame
  status=poc_getvar2d(filename,info.id,0,buf,mask,info);
  
  return status;
}


/*******************************************************************************
       TEMPLATES FOR NEW STUFF
*******************************************************************************/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_put_vara(int ncid,const cdfvar_t & info, size_t frame, const T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_put_vara(int,int,const size_t[],const size_t[],T*)
/*----------------------------------------------------------------------------*/
{
  int status,varid,iDI;//,,index of info.dim
  status=nc_inq_varid(ncid,info.name,&varid);

  size_t *start,*count;
  start=new size_t[info.ndim];
  count=new size_t[info.ndim];
  for(iDI=0;iDI<info.ndim;iDI++){
    if(isT(info.dim[iDI])){
      start[iDI]=frame;
      count[iDI]=1;
      }
    else{
      start[iDI]=0;
      count[iDI]=info.dim[iDI].length;
      }
    }
  
  status=poc_put_vara(ncid,varid,start,count,z);
  delete[]start;
  delete[]count;
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int poc_put_vara(const char *filename, const cdfvar_t &info, size_t frame, const T *z, int verbose=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for poc_put_vara(int,cdfvar_t,size_t,T*)
/*----------------------------------------------------------------------------*/
{
  int ncid,status;
  
  status=nc_open(filename,NC_WRITE,&ncid);
  if(status != NC_NOERR) goto error;
  
  status=poc_put_vara(ncid, info, frame, z);
  
  nc_close(ncid);
  
  error:
  if(status != NC_NOERR && verbose) nc_check_error(status,__LINE__,__FILE__,"error with %s on %s",info.name,filename);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int poc_gettime(T file, cdfvar_t variable, cdfgbl_t global, date_t *origine, double **time, int *nframes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Get times of frames of a given NetCDF file
See poc_gettime(const char*,cdfvar_t,cdfgbl_t,date_t*,double**,size_t*) for more help.
Declared because of size_t<->int conflict
*/
/*----------------------------------------------------------------------------*/
{
  size_t snf;
  int status;
  
  status=poc_gettime(file, variable, global, origine, time, &snf);
  *nframes=snf;
  
  return status;
}


extern int quoddy_loadr1(const char *name, int nndes, short  *buffer);
extern int quoddy_loadr1(const char *name, int nndes, float  *buffer);
extern int quoddy_loadr1(const char *name, int nndes, double *buffer);

// #error put poc_get_vara here

static int quoddy_loadr1(const char *name, int nndes, char *buffer) __attribute__((unused));
static int quoddy_loadr1(const char *name, int nndes, char *buffer) {
  TRAP_ERR_EXIT(ENOEXEC,"not coded for char\n");}

static int quoddy_loadr1(const char *name, int nndes, signed char *buffer) __attribute__((unused));
static int quoddy_loadr1(const char *name, int nndes, signed char *buffer) {
  TRAP_ERR_EXIT(ENOEXEC,"not coded for signed char\n");}


#endif
