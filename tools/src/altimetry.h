#if ALTIMETRY_H == 0
#define ALTIMETRY_H 1

typedef struct {
  int n;
  int *track,*segment;
  float *alias,*alongtrack,*harmonic;
  float *cut,*depth;
  complex<float> *zraw,*zlwf,*zhgf;
  complex<float> *smoothed;
  } xerror_t;

static float          data_rmask=1.e+10;
static complex<float> data_cmask=1.e+10;

extern int compute_coastal_error(char *wave, int kk, vector<mgr_t> mgr,int nmgr, xerror_t xerror, double *buble, int **list, int *card);
extern int compute_analysis_error(char *wave, vector<mgr_t> mgr,int nmgr,float *analysis_error);

#endif