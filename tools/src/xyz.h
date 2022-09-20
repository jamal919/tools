
#ifndef __XYZ
#define __XYZ

#define XYZ 0
#define YXZ 1

#define ROW    0
#define COLUMN 1

#define ASCII 0
#define SHAPE 1
#define CORA  2


#include "polygones.h"
#include "matrix.h"

extern  int reorder(double *x, int nvalues);

extern  int count_columns(char *s);
extern  int xyz_countdata (const char* filename);

extern int xyz_loadraw(const char *filename, string header, char *proj4_options, double * &x, double * &y, double * &z, double *mask, int & ndata, bool debug=false);
extern int xyz_load_shp(const char *filename, char *proj4_options, double * &x, double * &y, double * &z, double *mask, int & ndata);

extern int xyz_save(const char *filename, double *x, double *y, double **z, double mask, int ndata, int ncols);

extern int xyz_save(const char *filename, double *x, double *y, double *z, char *keep, int flag, int ndata);

extern int xyz_save(const char *filename, double *x, double *y, double *z, double mask, int ndata);
extern int xyz_save(const char *filename, double *x, double *y, double *z, double mask, int ndata, int header);
extern int xyz_save(const char *filename, double *x, double *y, short *z, short mask, int ndata);
extern int xyz_save(const char *filename, double *x, double *y, short *z, short mask, int ndata, int header);
extern int xyz_save(const char *filename, double *x, double *y, float *z, float mask, int ndata);

extern int xyz_save(const char *filename, grid_t grid, short *z,  short mask);
extern int xyz_save(const char *filename, grid_t grid, float *z,  float mask);
extern int xyz_save(const char *filename, grid_t grid, double *z, double mask);
extern int xyz_save(const char *filename, grid_t grid, double *z, double mask, int header);


extern  int xyz_readgrid_00(const char *,grid_t *, int);
extern  int xyz_readgrid_01(const char *,grid_t *, int, int);
extern  int xyz_readgrid_01_from_v1 (char *,grid_t *, int);


extern  int xyz_loadgrid(const char *, char *,grid_t *,int *);
extern  int xyz_loadmap(const char *,grid_t *);


extern int xyz_read00_r1(FILE *in, grid_t grid, int pos, int ncol, float *buf,float mask,int mode);

extern int xyz_read01_r1(FILE *in, grid_t grid, int pos, int ncol, float *buf,int mode);

extern int xyz_loadr1(const char *filename, grid_t & grid, int pos, float *buf, float *mask);
extern int xyz_loadr1_from_v1 (const char *, grid_t, float *, float *);

extern  int xyz_skipheader(FILE *in);

extern int xyz_PolygonSelection(vector<plg_t> & polygons, double* &x, double* &y, double* &z, int & ndata, frame_t & frame);
extern int xyz_PolygonSelection(plg_t *polygons, int npolygons, double* &x, double* &y, double* &z, int & ndata, frame_t & frame);
extern int xyz_PolygonSelection(const char     *polygons, double* &x, double* &y, double* &z, int & ndata, frame_t & frame);

extern int xyz_FrameSelection ( frame_t & frame, double* &x, double* &y, double* &z, int & ndata);

extern int xyz_CheckDuplicate_spherical(double *lon, double *lat, double *z, int ndata, char *keep, double resolution);
extern int xyz_CheckDuplicate_cartesian(double *lon, double *lat, double *z, int ndata, char *keep, double resolution);
extern int xyz_CheckDuplicate(double *t, double *p, double *z, int ndata, char *keep, double resolution, projPJ projection);

extern int xyz_decimate(frame_t frame, double *t, double *p, double *z,char *keep, int ndata, double resolution, double factor, bool cartesian, int verbose);

extern int xyz_reduce(double * &x, double * &y, double * &z, char *keep, int  &ndata);

#endif
