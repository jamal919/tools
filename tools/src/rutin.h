#include "tools-structures.h"

extern void gmerror(const char* error_text);
extern void gmerror(std::string error_text);
extern void bread_dg(char *addr, size_t size, size_t n, FILE *file);
extern void bwrite_dg(char *addr, size_t size, size_t n, FILE *file);
extern void free_imatrix(int **m,int nrl,int nrh,int ncl,int nch);
extern void free_ivector(int *v,int nl,int nh);
extern void free_smatrix(float **m,int nrl,int nrh,int ncl,int nch);
extern void free_dmatrix(double **m,int nrl,int nrh,int ncl,int nch);
extern void free_svector(float *v,int nl,int nh);
extern void free_dvector(double *v,int nl,int nh);
extern void gmplot(float yy[],unsigned short n,double x1,double dx,FILE *file);
extern int **imatrix(int nrl,int nrh,int ncl,int nch);
extern int *ivector(int nl,int nh);
extern float **smatrix(int nrl,int nrh,int ncl,int nch);
extern float *svector(int nl,int nh);
extern double **dmatrix(int nrl,int nrh,int ncl,int nch);
extern double *dvector(int nl,int nh);
