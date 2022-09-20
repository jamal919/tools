#if ARCHIVE_H == 0
#define ARCHIVE_H 1

#include "tools-structures.h"
#include "poc-time.h"
#include "map.h"

#include "poc-array.hpp"


// typedef struct {
//   int    type,nframe;
//   double scale[3];
//   double sampling,start;
//   date_t reference;
//   mesh_t mesh;
//   grid_t grid;
//   } meta_archive_t;

class meta_archive_t {
  private:
  public:
  int type; /* 1: byte with increment, 2: short, 4: float. See also clx_archiveinfo() */
#define FLAG_STATE 0
#define FLAG_FORCING 1
  int flag; /* 0: state vector, 1: forcing vector */
  int nframe;
  double scale[3];
  double sampling, start;
  date_t reference;
  mesh_t mesh;
  grid_t grid;

  meta_archive_t(){
    type = -1;
    flag = -1;
    nframe=-1;
    aset(scale,3,1.);
    sampling=0;
    start=-1;
//    reference=null_date;
    }

  void destroy(){
    mesh.destroy();
    }
};


typedef struct {
  float  **h1,**h2;
  double time1,time2;
  meta_archive_t info;
  int frame;
  char *file;
  int  activated;
  } memory_t;


/*--------------------------------------------------------------------------
 dummy_t  was previously define in the define of archive
---------------------------------------------------------------------------*/

class dummy_t {
  private:
  public:
  double time;
  meta_archive_t info;
  int frame;
  char *file;
  int  activated;
  
  dummy_t(){
    time=0;
    frame=0;
    file=0;
    activated=0;
    }
  } ;


extern int cdf_archiveread(const char *filename,const meta_archive_t & info, int frame, float **buffer, double *time);
extern int cdf_readexternal(double time, date_t reference, float *buffer[5], double *actual_time);
extern int cdf_initexternal(double time, double dT, state2D_t *set,  mesh_t mesh);
extern int cdf_archiveinfo(const char *filename, meta_archive_t *info);

extern void getnewdate(date_t reference, double t,date_t *actual);
extern char *sgetnewdate(date_t reference, double t);
extern int clx_initexternal(double time, double dT, state2D_t *set,  mesh_t mesh);
extern int clx_readexternal(double time, float *buffer[5], double *actual_time);
extern int clx_archiveread(FILE *file,const meta_archive_t & info, int step, float *buffer[3], double *time);

extern int clx_archivereadheader(FILE *file, meta_archive_t *info);
extern int clx_archivereadheader(const char *filename, meta_archive_t *info);

extern int clx_archivewriteheader(char *filename, meta_archive_t info, float **buffer);
extern int archive_createheader(char *path,meta_archive_t info,float **buffer);
extern int archive_skipheader(FILE *file, meta_archive_t info,double t, size_t size);

extern int archive_update2(char *filename, meta_archive_t info, float **buffer, double time);
extern int clx_archiveinfo(const char *filename, meta_archive_t *info);
extern int interpolate_external(basic_obsolete_t *basic, int n, double *out, float buffer[]);
extern int elements_statev1(meta_archive_t info, state2D_t *set, int nndes,mesh_t mesh);
extern int elements_statev2(meta_archive_t info,state2D_t *set,int nndes);
extern int instantaneous_external_state(double t, double dT, double step, state2D_t *set,int nndes);
extern void archive_setverbose(int flag);
extern int archive_read(char *filename, int fmt,const meta_archive_t & info, int step, float *buffer[3], float *mask, double *time);
extern int archive_freeinfo(meta_archive_t *info);
extern int archive_info(const char *filename, int fmt, meta_archive_t *info);

//from io.cpp
extern int section_loadr1(char *name, int nndes, float *buffer,int column);

extern int quoddy_loadr1(const char *name, int nndes, short  *buffer);

extern int quoddy_loadr1(const char *name, int nndes, float  *buffer, char *text);
extern int quoddy_loadr1(const char *name, int nndes, float  *buffer);
extern int quoddy_loadr1(const char *name, int nndes, double *buffer);
extern int quoddy_loadd1(const char *name, int nndes, double *buffer);

extern int quoddy_loadr2(const char *name, int nndes, float  *bufx, float  *bufy);
extern int quoddy_loadr2(const char *name, int nndes, double *bufx, double *bufy);

extern int quoddy_loadc1(const char *name, int nndes, complex<float> *buffer);
extern int quoddy_loadc1(const char *name, int nndes, complex<double> *buffer);

extern int quoddy_loadc2(char *name, int nndes, fcomplex *bufferx, fcomplex *buffery);
extern int load_s2r(char *filename, float *buf, int nnde,  int verbose);
extern int load_v2r(char *filename, float *bufx, float *bufy, int nnde,  int verbose);

extern int quoddy_saver1(const char *name, int nndes, float  *buffer, char **comment);
extern int quoddy_saver1(const char *name, int nndes, double *buffer, char **comment);

extern int quoddy_savec1(const char *name, int nndes, fcomplex *buffer, char **comment);
extern int quoddy_savec1(const char *name, int nndes, complex<double> *buffer, char **comment);

extern int quoddy_savec1(const char *name, int nndes, float *a, float *G, char **comment);

extern int quoddy_saver2(const char *name, int nndes, float  *bufferx, float  *buffery,char **comment);
extern int quoddy_saver2(const char *name, int nndes, double *bufferx, double *buffery,char **comment);

extern int quoddy_savec2(const char *name, int nndes, fcomplex *bufferx, fcomplex *buffery, char **comment);
extern int quoddy_savec2(const char *name, int nndes, dcomplex *bufferx, dcomplex *buffery, char **comment);

extern int quoddy_savec2(const char *name, int nndes, float *a_u, float *G_u, float *a_v, float *G_v, char **comment);
extern int quoddy_savec2(const char *name, int nndes, double *a_u, double *G_u, double *a_v, double *G_v, char **comment);

extern int fe_framemapr1(char *input, meta_archive_t info, grid_t grid, float mask, int *elts, int step, float *buf, double *tag);
extern int fe_framemapr2(char *input, meta_archive_t info, grid_t grid, float mask, int *elts, int step, float *bufx, float *bufy, double *tag);

#define TUGO_UG_BINARY_FORMAT 0
#define TUGO_UG_NETCDF_FORMAT 1
#define TUGO_SG_NETCDF_FORMAT 2

#endif
