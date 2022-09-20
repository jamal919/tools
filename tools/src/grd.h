
#include "tools-structures.h"

#ifndef GRD_H
#define GRD_H
extern int grd_loadgrid (const char* filename,grid_t *grid);
extern int grd_savegrid (const char* filename, grid_t grid);
extern void grd_getdimension(char *filename,int* nv,int* ni,int* nj,int* nk,int* nd,int* nt,int* status);
extern int grd_checkfile (char *filename);

extern int grd_loadr1   (const char* filename, grid_t grid, float *buf, float *mask);
extern int grd_loadr1   (const char* filename, grid_t grid, short *buf, short *mask);

extern int grd_loads1   (const char* filename, grid_t grid, int n, short *buf, short *mask);

// extern void grdf_loadr1 (const char* filename, int *k,int *d,int *t,grid_t *grid, int *n, float *buffer, float *mask, int *status);

// extern float *grd_mirror_r_obsolete (grid_t grid, int n, float *buf, float mask);
extern int grd_mirror_r (grid_t grid, int n, float *buffer, float mask);
extern int grd_mirror_r (grid_t grid, int n, short *buffer, short mask);
extern int grd_mirror (grid_t grid, int n, complex<float> *buffer, complex<float> mask);

extern int  grd_extract (const char* filename, grid_t grid, int n, float *buf, float *mask);
//extern void grdf_extractr1 (const char* filename, grid_t *grid, int *n, float *buffer, int *shift, int *status);
extern int  grd_save      (const char* filename, grid_t grid, int n, float *buf, float mask);
// extern void grdf_saver1   (const char* filename, int *k,int *d,int *t, int *tmax,grid_t *grid, int *n, float *buffer, float *mask, int *status);
extern int  grd_save      (const char* filename, grid_t grid, int n, short *buf, short);
// extern void grdf_saves1   (const char* filename,grid_t *grid, int *n, short *buffer, int *status);
// extern void pivote_buffer (double *buffer,grid_t *grid,double mask);

extern int grd_DuplicateTag(const char* in, const char* out);
#endif
