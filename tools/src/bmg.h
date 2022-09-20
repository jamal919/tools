#include "tools-structures.h"
#include "map.h"

extern int filesize (char *cmd);
extern void bmg_printheader (bmgheader_t header);
extern int bmg_readheader (FILE *in,bmgheader_t *header);
extern int bmg_writeheader (FILE *in,bmgheader_t header);
extern int bmg_createheader (char* file,int kmax,int dmax,int tmax,grid_t grid,float *levels,float mask);
extern int bmg_read_r (FILE *in, bmgheader_t header, int k,int d,int t, float *buf, float *time);
extern void bmgf_read_r (FILE *in, bmgheader_t *header, int *k,int *d,int *t, float *buf, float *status);
extern int bmg_read_d (FILE *in, bmgheader_t header, int k,int d,int t, double *buf, float *time);
extern int bmg_write_r (FILE *in, bmgheader_t header, int k,int d,int t, float *buf, float time);
extern int bmg_loadgrid (const char* filename,grid_t *grid);
extern int bmg_getinfo (char* filename, bmgheader_t *header);
extern void bmg_h2g (bmgheader_t header, grid_t *grid);
extern void bmg_checkfile (char *filename, int *status);

extern int bmg_loadr1 (const char* file,int k,int d,int t,grid_t grid,float *buf,float *mask,float *time);
extern int bmg_loadr1 (const char* file,int k,int d,int t,grid_t grid,short *buf,short *mask,float *time);

extern int bmgf_loadr1 (char* file, int *k,int *d,int *t,grid_t *grid,float *buf,float *mask,float *time);
extern int bmg_save_r (char* file, int t, grid_t grid, float *buf, float time);
extern int bmgf_save_r(char* file, int *t,grid_t *grid,float *buf,float *time);
extern int bmg_saver1 (char* file, int k, int d, int t, grid_t grid, float *buf, float time,float mask);
extern int bmgf_saver1(char* file, int *k,int *d,int *t,grid_t *grid,float *buf,float *mask,float *time);
extern int bmg_saver2 (char* file, int k, int d, int t, grid_t grid, float *bufx, float *bufy, float time, float mask);
extern int bmg_loadc1(const char* file,int k,int d,int t,grid_t *grid,fcomplex **buffer,fcomplex *mask,float *time);
extern int bmgf_loadc1(char* file,int *k,int *d,int *t,grid_t *grid,fcomplex *buf,fcomplex *mask,float *time);
extern int bmg_savec1 (char* filename,int k,int t,int d,grid_t grid,fcomplex *buf,float time,fcomplex mask);
extern int bmg_savec2 (char* filename,int k,int t,int d,grid_t grid,fcomplex *bufx,fcomplex *bufy,float time,fcomplex mask);
