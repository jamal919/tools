
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief structured grid function declarations
*/
/*----------------------------------------------------------------------------*/

#if MAP_H == 0
#define MAP_H 1

#include "poc-time-classes.h"

#include "maths.h"
#include "geo.h"
#include "statistic.h"

#include "map-classes.h"


extern frame_t limits(const grid_t *grids,int ngrids);
extern frame_t limits(const grid_t & grid);

extern float    map_recale(const grid_t & grid,float x);
extern double   map_recale(const grid_t & grid,double x);
extern void map_recale_all_x(grid_t *grid);

extern int *map_get_perimeter_indexes(const grid_t & grid,int *n,const double *field=0,double mask=NAN,int detectMaskedLine=0,int verbose=0);

extern grid_t   map_getgrid2d(const grid_t & grid3d);
extern grid_t   map_getgrid3d(const grid_t & grid2d) ;

extern int mapc_checkconnexity(const grid_t & grid);

extern statistic_c_t cget_geostatistics(const grid_t & grid,complex<float> *h, complex<float> mask);

extern void map_allocate_x_y(grid_t *grid,int modeH=-1,int nx=-1,int ny=-1);

extern int map_set_dxdy(grid_t *grid);
extern int map_completegridaxis_2(grid_t *grid, int transposed=0);
extern int map_completegridaxis(grid_t *grid,int mode=2,int verbose=0);
extern int map3d_completegridaxis(grid_t *grid);

//extern int map_loadfield_cdf(const char *input,const char *varname, const char *xname,const char *yname, const grid_t & grid, float * &buffer, float & mask);
extern int map_loadfield_cdf(const char *input,const char *varname, grid_t & grid, float * &buffer, float & mask);

extern int map_set2Dgrid(grid_t *grid, double xmin, double ymin, double xmax, double ymax, double dx=0., double dy=0.);
extern int map_minmax(grid_t *grid);
extern int map_set2Dgrid(grid_t *grid, frame_t frame, double dx=0., double dy=0.);
extern int map_set2Dgrid(grid_t *grid, frame_t frame, int nx, int ny);
extern int map_set2Dgrid(grid_t *grid, frame_t frame, int n);

extern grid_t get_zonegrid(const char *zone, bool *identified=0);
extern bool read_zone_arg(const char *keyword,const char *arg,frame_t *prescribed,double *dx,double *dy);
extern void print_zone_arg_help();
extern int apply_zone_arg(grid_t *grid,const frame_t &prescribed,double dx,double dy);

extern int map_cartesian_resolution(const grid_t & grid, double **dx, double **dy, double **dz=NULL);
extern int map_spherical_resolution(const grid_t & grid, double **dx, double **dy, double **dz=NULL);

extern double map_resolution(const grid_t & grid, double t, double p, double *d, double *dy);
extern double map_resolution(const grid_t & grid, int k, int l, double & dx, double & dy, int mode);

extern float map_getf_xyt(const grid_t & grid, float *buf1, float *buf2, double t1, double t2, float mask,float x, float y, double t);

extern int map_coeffcients02(const grid_t & grid, double x, double y, double *area, int *node, int *count);
extern int map_interpolation3d(const grid_t & grid,double x, double y, double z, float *buffer);
extern int map_index00(const grid_t & grid, double x, double y, int *ktrue, int *ltrue);
extern int map_index02(const grid_t & grid, double x, double y, int *ktrue, int *ltrue);

extern bool map_check_circular(const grid_t & grid);
extern void set_grid_list(grid_t *grid,int wrap=0,int verbose=0);

extern int ind1D(const int    *t,int n,int    time,int *m);
extern int ind1D(const float  *t,int n,float  time,int *m);
extern int ind1D(const double *t,int n,double time,int *m);
extern int ind1D(const vector<double> & t,int n,double time,int *m);
extern void index_interpolation(const grid_t &grid,double lon,double lat,int64_t *m,const float *buf,const float & mask,float *z,int verbose=0);
extern void index_interpolation(const grid_t &grid,double lon,double lat,int64_t *m,const double *buf,const double & mask,double *z,int verbose=0);
extern void index_interpolation(const grid_t &grid,double lon,double lat,int64_t *m,const double *buf,const double & mask,complex<double> *z,int verbose=0);
extern void index_interpolation(const grid_t &grid,double lon,double lat,int64_t *m,const complex<double> *buf,const complex<double> & mask,complex<double> *z,int verbose=0);
extern void index_interpolation(const grid_t &grid,double lon,double lat,int64_t *m,const complex<float> *buf,const complex<float> & mask,complex<float> *z,int verbose=0);

extern int map_index(const grid_t & grid, double x, double y, int *ktrue, int *ltrue);

extern int map_interpolation(const grid_t & grid,signed char *dum,signed char mask,double x,double y, float *z);
extern int map_interpolation(const grid_t & grid,short  *dum, short mask,double x,double y, float *z);
extern int map_interpolation(const grid_t & grid,float  *dum, float mask,double x,double y, float *z);
extern int map_interpolation(const grid_t & grid,double *dum,double mask,double x,double y,double *z);
extern int map_interpolation(const grid_t & grid,complex<float> *buf,complex<float> mask,double x,double y,complex<float> *z);
extern int map_interpolation(const grid_t & grid,double *dum,double mask,double x,double y,complex<double> *z);
extern int map_interpolation(const grid_t & grid,complex<double> *buf,complex<double> mask,double x,double y,complex<double> *z);

extern int map_average (const grid_t & grid, float *buffer, float mask, grid_t target, int m0, int n0, float *z);
extern int map_average (const grid_t & grid, double *buffer, double mask, grid_t target, int m0, int n0, double *z);
extern int map_average (const grid_t & grid, complex<float> *buffer, complex<float> mask, grid_t target, int m0, int n0, complex<float> *z);
extern int map_average (const grid_t & grid, complex<double> *buffer, complex<double> mask, grid_t target, int m0, int n0, complex<double> *z);

extern void distance_to_nearest_unmasked(const grid_t &grid,const double *buffer,double mask,double *distance,double *value=0);
extern void distance_to_nearest_unmasked(const grid_t &grid,const complex<double> *buffer,complex<double> mask,double *distance,complex<double> *value=0);
extern void distance_to_nearest_unmasked(const grid_t &grid,const float *buffer,float mask,double *distance,float *value=0);
extern void distance_to_nearest_unmasked(const grid_t &grid,const char *buffer,char mask,double *distance,char *value=0);

extern int map_export(const grid_t & g_source, float *b_source, float m_source,grid_t g_target, float *b_target, float m_target,int mode=-1);
extern int map_export(const grid_t & g_source, complex<float>  *b_source, complex<float>  m_source,grid_t g_target, complex<float>  *b_target, complex<float>  m_target,int mode=-1);
extern int map_export(const grid_t & g_source, complex<double>  *b_source, complex<double>  m_source,grid_t g_target, complex<double>  *b_target, complex<double>  m_target,int mode=-1);

extern int map_remap(grid_t & g_source, complex<float> *b_source, complex<float> m_source, grid_t *g_target, complex<float> **b_target, complex<float> *m_target,int incr,int mode);
extern int map_remap (grid_t & g_source, float *b_source, float m_source, grid_t *g_target, float **b_target, float *m_target, int incr,int mode);

extern int map_gradient(const grid_t & grid,int n, float *dum, float mask, int ref, float *dumx, float *dumy);
extern int map_gradient(const grid_t & grid,int n, complex<float> *dum, complex<float> mask, int ref, complex<float> *dumx, complex<float> *dumy);

// extern int map_curve00(const grid_t & grid,int n, float *dum, float mask, int ref, float *dumx, float *dumy);
// extern int map_curve02(const grid_t & grid,int n, float *dum, float mask, int ref, float *dumx, float *dumy);

extern int map_curve(const grid_t & grid, float *dum, float mask, int ref, float *dumx, float *dumy);
extern int map_curve(const grid_t & grid, double *dum, double mask, int ref, double *dumx, double *dumy);
extern int map_curve(const grid_t & grid, complex<double> *dum, complex<double> mask, int ref, complex<double> *dumx, complex<double> *dumy);

extern int map_laplacian(const grid_t &grid, complex<double> *dum, complex<double> mask, int ref, complex<double> * & laplacian);
extern int map_divergence(const grid_t &grid, double *dumx, double *dumy, double mask, int ref, double *divergence);

extern void map_printgrid (const grid_t & grid);
extern void mapc_printgrid3d (const grid_t & grid);
extern int wheader(FILE *out,char **comment, grid_t grid, float spec, int n, float *average);
extern int rheader( FILE *out, grid_t *grid, float *spec, int *M, int *n, float *average);
extern int mapc_index(const grid_t & grid, double x, double y, int *ktrue, int *ltrue);

extern int map_interpolate1D(const float *h,const float *t,  float mask,  int n, float time,  float *z,  int extrapolate=0, int *k0=0);
extern int map_interpolate1D(const float *h,const double *t, float mask,  int n, double time, float *z,  int extrapolate=0, int *k0=0);
extern int map_interpolate1D(const double *h,const double *t,double mask, int n, double time, double *z, int extrapolate=0, int *k0=0);
extern int map_interpolate1D(const double *h,const double *t,double mask,int n,double time,complex<double> *z,int extrapolate=0, int *k0=0);
extern int map_interpolate1D(const complex<double> *h,const double *t,complex<double> mask,int n,double time,complex<double> *z,int extrapolate=0, int *k0=0);
extern int map_interpolate1D(const double *h,const double *t, double mask, int n, double time, double *z, double tmax, int extrapolate=0, int *k0=0);
extern int map_interpolate1D(const double *h,const double *t, double mask, int n, double time, short unsigned int *z, double tmax, int extrapolate=0, int *k0=0);
extern int map_interpolate1D(const vector3_t *h,const double *t, int n, double time, vector3_t *z, int *k0);

extern int map_interpolate1D(const std::complex<float> *h,const double *t, std::complex<float> mask, int n, double time, std::complex<float> *z);
extern int map_interpolate1D(const vector<float> & h,const vector<double> & t, float mask, int n, double time, float *z, int extrapolate=0, int *k0=0);


extern int map_nearestvalue00(const grid_t & grid,int n, float *buf, float mask, double x, double y, float *z, int maxdepth=10);
extern int map_nearestvalue00(const grid_t & grid, int n, unsigned short *buf, unsigned short mask, double x, double y, unsigned short *z, int maxdepth=10);

extern int map_ClosestVertex(const grid_t &grid, double x, double y, int & m0, double & distance, bool debug);

extern int check_vertical_direction(const double *z, double zmask, size_t HSize, int nz, bool *increasing, double *factor);
extern int check_vertical_direction(const grid_t & grid, bool *increasing, double *factor);
extern int map_sliceH(const grid_t & grid, double level, float *buffer, float mask, float *slice);
extern int map_sliceH(const grid_t & grid, double level, double *buffer, double mask, double *slice);
extern int map_sliceH(const grid_t & grid, double level, double *buffer, double mask, complex<double> *slice);
extern int map_sliceH(const grid_t & grid, double level, complex<double> *buffer, complex<double> mask, complex<double> *slice);
extern int map_sliceV(const grid_t & grid, double *x, double *y, int nloc, float *buffer, float mask,grid_t *slice_grid, float **slice, int verbose=1);
extern int map_sliceV(const grid_t & grid, double *x, double *y, int nloc, double *buffer, double mask,grid_t *slice_grid, double **slice, int verbose=1);
extern int map_sliceV(const grid_t & grid, double *x, double *y, int nloc, double *buffer, double mask,grid_t *slice_grid, complex<double> **slice, int verbose=1);
extern int map_sliceV(const grid_t & grid, double *x, double *y, int nloc, complex<double> *buffer, complex<double> mask,grid_t *slice_grid, complex<double> **slice, int verbose=1);
extern int map_profile01(const grid_t & grid, double x, double y, float *buffer, float mask,grid_t *slice_grid, float **profile_x, float **profile_y, int *nloc);
extern int map_chklevels(const grid_t & grid);
extern int map_diffusion01(const grid_t & grid, float *buffer, float mask);
extern int map_profile02(const grid_t & grid, double x, double y, float *buffer, float mask,float **profile_x, float **profile_y, int *nloc);

extern int map_sliceV3D(const grid_t & grid, double *x, double *y, int nloc, float *buffer, float mask, grid_t *slice_grid, float **slice, float **bottom,int xoption,int yoption);
//extern int map_sliceV3D(const grid_t & grid, double *x, double *y, int nloc, double *buffer, double mask, grid_t *slice_grid, double **slice, double **bottom,int xoption,int yoption);
extern int map_sliceV3Dvector(const grid_t & grid, double *x, double *y, int nloc, float *bufferx, float *buffery, float mask,grid_t *slice_grid, float **slice, float **bottom,int xoption, int yoption, int voption);

extern double map_grid_x (const grid_t & grid, int k, int l);
extern double map_grid_y (const grid_t & grid, int k, int l);
extern int map_integraleR(const grid_t & grid, float *buf, float mask, double *sum,double *area);
extern int map_integraleC(const grid_t & grid, complex<float> *buf, complex<float> mask, complex<double> *sum,complex<double> *area);

//extern int copy_grid_to_grid1d(const grid_t & grid3d,grid_t *grid1d);

extern int map_loadfield_obsolete(const char *filename,grid_t *grid, float **buffer, float *mask);
extern int map_loadfield(const char *filename,const char *varname, grid_t *grid, float **buffer, float *mask);

extern int map_loadfield(const char *filename,int format, const char *varname, grid_t *grid, float **buffer, float *mask);
extern int map_loadfield(const char *filename,int format, const char *varname, grid_t *grid, short **buffer, short *mask);

extern int map_loadfield(const char *filename,int format, const char *varname, const char *xname, const char *yname, const char *proj, grid_t *grid, float **buffer, float *mask);
extern int map_loadfield(const char *filename,int format, const char *varname, const char *xname, const char *yname, const char *proj, grid_t *grid, short **buffer, short *mask);

extern  int map_extendfield(grid_t *grid, complex<float> **buffer, int *modified);
extern  int map_extendfield(grid_t *grid, complex<float> **buffer, int nbuffers, int *modified);

extern  int map_persistence(const grid_t & grid, signed char *buf, signed char mask, double radius);
extern  int map_persistence(const grid_t & grid, float *buf, float mask, double radius);
extern  int map_persistence(const grid_t & grid, double *buf, double mask, double radius);
extern  int map_persistence(const grid_t & grid, complex<float> *buf, complex<float> mask, double radius);

extern  int map_smooth_local(const grid_t & grid, float *buf, float mask, float scale, int k, int l, float *out);
extern  int map_smooth_local(const grid_t & grid, complex<float> *buf, complex<float> mask, float scale, int k, int l, complex<float> *out);

extern  int map_smooth(int mode, const grid_t & grid, const float *buf, const float mask, const float scale, float *out);
extern  int map_smooth(int mode, const grid_t & grid, const double *buf, const double mask, const float scale, double *out);
extern  int map_smooth(int mode, const grid_t & grid, const complex<float> *buf, const complex<float> mask, const float scale, complex<float> *out);

extern  int map_smooth(const grid_t & grid, const float  *buf, const float  mask, const float *scale, float  *out);
extern  int map_smooth(const grid_t & grid, const double *buf, const double mask, const float *scale, double *out);

extern  int map_smooth_latThenLong(const grid_t & grid, const complex<float> *buf, const complex<float> mask, const float scale, complex<float> *out);
extern  int map_smooth_latThenLong(const grid_t & grid, const float *buf, const float mask, const float scale, float *out);

extern  int map_diffusion_forward(const grid_t & grid, const complex<float> *buf, const complex<float> mask, const float scale, complex<float> *out, int niteration);
extern  int map_diffusion_forward(const grid_t & grid, const float *buf, const float mask, const float scale, float *out, int niteration);

extern  int map_LengthScale(const grid_t & grid, complex<float> *tides, complex<float> mask, float* & lamda1, float* & lamda2, float* & lamda3);

extern  int map_projection_forward( grid_t & grid, const char *proj4);
extern  int map_projection_backward(grid_t & grid, const char *proj4);

extern  grid_t map_get_spherical(projPJ projection,grid_t sgrid);
extern  grid_t map_get_cartesian(projPJ projection,grid_t sgrid);

extern  grid_t map_get_spherical(geo_t projection,grid_t sgrid);
extern  grid_t map_get_cartesian(geo_t projection,grid_t sgrid);

extern  grid_t map_Rgrid(frame_t frame, double dx, double dy, int mode);
extern  grid_t map_Rgrid(frame_t frame, size_t nx, size_t ny, int mode);

extern  grid_t map_duplicategrid(grid_t zgrid);

extern  float *map_duplicate2D(grid_t zgrid, float *buffer, grid_t extended);
extern  float *map_duplicate3D(grid_t zgrid, float *buffer, grid_t extended);

extern  int parse_GridOptions(const string & options, metagrid_t & meta);

extern grid_t map_f2vgrid(const grid_t & fgrid);
extern grid_t map_f2ugrid(const grid_t & fgrid);

extern  float *map_interpolate_v2z(grid_t zgrid, grid_t vgrid, float *v_field, float mask);
extern  float *map_extrapolate_z2v(grid_t zgrid, grid_t vgrid, float *z_field, float mask);


extern  grid_t map_vgrid(const grid_t & zgrid, float *ztopo=0, float topomask=0, const char *output=0);

extern  grid_t map_z2ugrid(const grid_t & fgrid);
extern  grid_t map_z2vgrid(const grid_t & fgrid);
extern  grid_t map_f2zgrid(const grid_t & fgrid);

extern resize_t map_CheckZgrid(grid_t & zgrid, char *zlandmask);

extern  int map_resize_grid(grid_t & zgrid, resize_t resize);
extern  int map_resize(grid_t & zgrid, const resize_t & resize, char   * & buffer);
extern  int map_resize(grid_t & zgrid, const resize_t & resize, float  * & buffer);
extern  int map_resize(grid_t & zgrid, const resize_t & resize, double * & buffer);
extern  int map_get_format(const char *format);

#endif
