#include "tools-structures.h"

// function from mss-comp.v2.cpp

extern void trace( int nt, double *tmtna, double lambda0, double *lon, double *lat);
extern void interpTrace( char tTrace[3], double lonRef, double latRef, int nt, double *dt, double *lon, double *lat );
extern int compute_nominal_track(grid_t *backbone,grid_t *nominal_track);
extern int set_topographic(double *Cm, int dimM, grid_t mssgrid, float *mss_cls);
extern int set_gaussian(double *Cm, int dimM, grid_t mssgrid);
extern int set_byderivative(double *Cm, int dimM, grid_t mssgrid);
extern int compute(grid_t mssgrid,serie_t TPdata,float *mss_cls,double *error,float *innovation,float *analysis,float *G);
int compute_annual(float* h, double *t, float mask, int count);
grid_t get_trackgridsioux(double *lon,double *lat, serie_t data,grid_t *nominal_track,char sens[3]);
void get_extremities(char *sat,int k,serie_t TPdata,double *lon,double *lat);

//function from mss-io.v2.cpp


int mss_createfile(char *filename,size_t x_len, size_t y_len, grid_t grid);
int writefile(char *filename, float *mss, float *cls, float *ino);
int savenodes(char *filename, serie_t metadata);
serie_t load_metadata_ref(char *input, int point);
int decale_metadata(data_t *data,int *keep,int nmes);
serie_t load_metadata_raw(char *input);
void write_list_header(FILE *out,int n);
void save_meantracks_CLS(serie_t TPdata,char *track,float *mss,plg_t *polygones,int npolygones);
void save_meantracks_CTOH(serie_t TPdata,char *alti_file, int pair,float *mss,grid_t *nominal_track);




   //Variables globales pour les deux functions suivantes*/
   // a ameliorer car pas de variable globale dans un code CPP

