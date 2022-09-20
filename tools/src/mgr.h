/** \file
tide gauges stuff
*/

#ifndef MGR_H
#define MGR_H 1


#include "tools-structures.h"
#include "tides.h"
#include "polygones.h"


#define MGR_FORMAT_CODE_LEGOS_ASCII    1
#define MGR_FORMAT_CODE_LEGOS_NETCDF   2
#define MGR_FORMAT_CODE_LEGOS_OBS      3
#define MGR_FORMAT_CODE_DART           4
#define MGR_FORMAT_CODE_RAY            5
#define MGR_FORMAT_CODE_ATG_DATABASE   6
#define MGR_FORMAT_CODE_SHOM           7
#define MGR_FORMAT_CODE_BIO            8

#define MGR_FORMAT_NAME_LEGOS_ASCII    "LEGOS-ASCII"
#define MGR_FORMAT_NAME_LEGOS_NETCDF   "LEGOS-NETCDF"
#define MGR_FORMAT_NAME_LEGOS_OBS      "LEGOS-OBS"
#define MGR_FORMAT_NAME_DART           "DART"
#define MGR_FORMAT_NAME_RAY            "RAY"
#define MGR_FORMAT_NAME_ATG_DATABASE   "ATG_DATABASE"
#define MGR_FORMAT_NAME_SHOM           "SHOM"
#define MGR_FORMAT_NAME_SHOM           "BIO"


//from mgr.cpp
extern int mgr_loadobs(const char *filename, char *wave, int *nmgr ,vector<mgr_t> & set, double ***covariance, int ***list, int **card);

extern int mgr_save_obs(const char *filename, vector<mgr_t> &, const char *wave);
extern int mgr_save_obs4assim(const char *filename, char *poly, char *wave, int nmgr ,vector<mgr_t> , float **covariance, int **list, int *card);

extern int mgr_Constants2Legend_00(const char *filename, char *wave, int nmgr ,vector<mgr_t> , float **errors, float **covariance, int *track, int **list, int *card);
extern int mgr_Constants2legend_01(const char *filename, char *wave, vector<mgr_t>);

extern int mgr_Timeserie2Legend(const char *filename, double lon, double lat, double *time, int nframes, double **values, int nvalues);

extern map<string,int> MgrHarmonic_formats;
extern void print_mgrh_formats(void);
extern int init_mgrh_formats(void);

extern int mgr_load(const char *filename, vector<mgr_t> & mgr_serie, bool reset=true);
extern int mgr_load(const char *filename, vector<mgr_t> & mgr_serie, int format, bool reset=true);
extern int mgr_load(const char *filename, vector<mooring_t> & moorings, vector<spectrum_t> & spectrum, vector<hconstant_t> & constants, int format);


extern int mgr_load_ascii(const char *filename, vector<mgr_t> & mgr_serie);
extern int mgr_loadRAY(const char *filename, vector<mgr_t> & mgr_serie);
extern int mgr_loadATG_DATABASE(const char *filename, vector<mgr_t> & mgr_serie);

extern int mgr_save(const char *filename, vector<mgr_t> & mgr_serie, float);

extern int mgr_save_ascii(const char *filename, const vector<mgr_t> & mgr_serie);

extern int mgr_save2D(const char *filename, int nmgr ,vector<mgr_t> *mgr_serie, float);

extern int mgr_create(int nmgr ,vector<mgr_t> & mgr_serie, const spectrum_t& s);
extern int mgr_create(int nmgr, vector<mgr_t> & mgr, const spectrum_t& s, double *lon, double *lat, double *depth);

extern int mgr_save_netcdf(const char *filename ,vector<mgr_t> _serie);
extern int mgr_load_netcdf(const char *filename, vector<mgr_t> & mgr_serie);

extern int mgr_read_ray(char *filename, int *nmgr ,vector<mgr_t> & mgr_serie);

extern int mgr_save_BIO(const char *filename , vector<mgr_t> & mgr_serie);

extern int mgr_save(const char *filename, vector<mgr_t> mgr, const string & format, const char *wave=0);


extern int *mgr_select(const char *poly ,vector<mgr_t> mgr);
extern int *mgr_select(vector<plg_t> & polygons, vector<mgr_t> & mgr);
extern int mgr_select(vector<mgr_t> , int *list, unsigned short *, unsigned short);

extern  int mgr_exclude(vector<mgr_t> & mgr, const char *poly, int position=PLG_POINT_EXTERIOR);
extern  int mgr_exclude(vector<mgr_t> & mgr, vector<plg_t> & polygons, int position=PLG_POINT_EXTERIOR);

extern  int mgr_concat(vector<mgr_t> & mgr, vector<mgr_t> additionals);

typedef struct {
  int id;
  int node[4];
  int track[2];
  double t,p,d;
//  double a[2],G[2];
//  std::complex<float> c[2];
//  double a[2],G[2];
  std::complex<float> *c[2];
  std::complex<float> *error;
  } xover_t;

typedef struct {
  int n;
  xover_t *xover;
  } xoverbase_t;

extern int xover_id(xoverbase_t base,xover_t *xover);
extern xoverbase_t read_xover(char *input);

extern int *mgr_order(vector<mgr_t> ,int nmgr);
extern int *mgr_order(vector<mgr_t> ,int nmgr,const char *criterion, int direction);
extern int mgr_WaveOrder(vector<mgr_t> _serie, int nmgr, string orderCrit);

extern int mgr_compact(vector<mgr_t> & mgr, int& nmgr);
extern string shellGetDate();
extern void createFile(const char *filename, int nWaves, int nRecords);
extern void putVariables(const char *filename, int nmgr, vector<mgr_t> _serie);

//from mgr-woce.cpp
extern int hawaii_decode_position_p(char *line, mooring_t *local, int *year=0);
extern date_t hawaii_decode_data(FILE *in, char *line);
extern int mgr_loadGLOSS_local(const char *filename, date_t first, date_t last, mooring_t *mooring, double **elevation, double **time, double *mask, int *n, const char *outfile);
extern int hawaii_getlength(const char *filename,date_t *start,date_t *finish,mooring_t *mooring);
extern int read_hawaii_positions(const char* input, pressure_station** sample);
extern int read_hawaii(FILE *hawaii,double *lat,double *lon,double *time,double *sealevel,date_t firstdate,int n);
extern void hawaii_write_header(const char* strcode, int year, double lat, double lon, FILE* hawaii);
extern int hawaii_write_data(const char* strcode, double time, double* data, double mask, FILE* hawaii, double lat, double lon, int refyear);
extern void mgr_savehawaii(FILE* hawaii, double lat, double lon, double* residual, double* time, int numstation, const char* stationname, date_t firstdate, int nmes, double t1);
extern int mgr_savehawaii(const char *filename,double lat,double lon,double *h,double *t,int code,char *name,date_t firstdate,int nvalues,double t1);

extern int mgr_loadGLOSS(const char *filename, tseries_t *serie, mooring_t *mooring);
extern int mgr_loadHAWAII_NC(const char *filename, tseries_t *serie, mooring_t *mooring);

//from mgr-moorings-stations.cpp
extern int mgr_read_moorings(const char *, mooring_t **);
extern int mgr_save_moorings(const char *, int , mooring_t *);

extern int mgr_stations_2_moorings(int nst, pressure_station *, mooring_t **);
extern int mgr_moorings_2_stations(int nst, mooring_t *, pressure_station **);

extern int mgr_read_stations(const char *, pressure_station **);
extern int mgr_save_stations(const char *, int , pressure_station *);

#endif
