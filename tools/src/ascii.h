#include "tools-structures.h"
#include "map.h"
extern void ascii_printheader (bmgheader_t header);
extern int ascii_readheader (FILE *in,bmgheader_t *header);
extern int ascii_writeheader (FILE *in,bmgheader_t header);

extern int  ascii_createheader (const char* file,int kmax,int dmax,int tmax,grid_t grid,float *levels,float mask);
extern void ascii_loadgrid     (const char* filename,grid_t *grid, int *status);
extern int  ascii_getinfo      (const char* filename, bmgheader_t *header);
extern void ascii_checkfile    (const char *filename, int *status);

extern int ascii_load (const char* filename, grid_t *grid,float **buf,float *mask);
extern int ascii_load (const char* filename, grid_t *grid,short **buf,short *mask);

extern int ascii_load_GIS (const char* filename, grid_t *grid, const char *proj, float **buf, float *mask);
extern int ascii_load_GIS (const char* filename, grid_t *grid, const char *proj, short **buf, short *mask);

extern int ascii_load_SLIM (const char* filename, grid_t *grid, const char *proj, float **buf, float *mask);
extern int ascii_load_SLIM (const char* filename, grid_t *grid, const char *proj, short **buf, short *mask);

extern int ascii_save_r (FILE *out, int t, grid_t grid, float *buf, float mask, double time);
extern int ascii_saver1 (const char* filename, grid_t & grid, float *buf, float mask, const char *format, int step);

extern int ascii_save_GIS(const char* filename, grid_t & grid, float *buf, float mask, const char *format, int step);

extern int ascii_loadc2 (const char* filename,grid_t *grid,complex<float> **buf1,complex<float> **buf2,complex<float> *mask);
extern int ascii_savec1 (const char* filename,grid_t grid,complex<float> *buf,complex<float> mask);
extern int ascii_savec2 (const char* filename,grid_t grid,complex<float> *bufx,complex<float> *bufy, complex<float> mask);

extern int ascii_loadc1_got (const char* filename,grid_t *grid,complex<float> **buffer,complex<float> *mask);
extern int ascii_loadc1_got_complex (const char* filename,grid_t *grid,complex<float> **buffer,complex<float> *mask);

extern int shomcst_loadgrid (const char* filename,grid_t *grid, int mode);
extern int shomcst_read (const char *filename, grid_t grid, int pos, float *buf, float *mask);
extern int shomcst_read (const char *filename, grid_t grid, int pos, complex<float> *buf, complex<float> *mask);

extern int shomcst_loadc1 (const char *filename, grid_t *grid, int pos, complex<float> ***buf, complex<float> *mask, spectrum_t *s);
