
#include "tools-structures.h"
#include <vector>

#ifndef XTRACK_IO
#define XTRACK_IO

#define X_TRACK_RAW_NUMVAL 31

#define X_TRACK_REF_TIME 1
#define X_TRACK_REF_SSHA 2
#define X_TRACK_REF_DAC_GLOBAL 3
#define X_TRACK_REF_DAC_GLOBAL_S1 4
#define X_TRACK_REF_DAC_GLOBAL_S2 5
#define X_TRACK_REF_DAC_REGIONAL 6
#define X_TRACK_REF_DAC_REGIONAL_S1 7
#define X_TRACK_REF_DAC_REGIONAL_S2 8
#define X_TRACK_REF_DAC_IB 9
#define X_TRACK_REF_DAC_IB_S1 10
#define X_TRACK_REF_DAC_IB_S2 11
#define X_TRACK_REF_DAC_LSA 12
#define X_TRACK_REF_DAC_SET 13
#define X_TRACK_REF_TIDE_REGIONAL 14
#define X_TRACK_REF_TIDE_REGIONAL_S1 15
#define X_TRACK_REF_TIDE_REGIONAL_S2 16
#define X_TRACK_REF_TIDE_GOT4_7 17
#define X_TRACK_REF_TIDE_GOT4_7_S1 18
#define X_TRACK_REF_TIDE_GOT4_7_S2 19
#define X_TRACK_REF_TIDE_FES2004 20
#define X_TRACK_REF_TIDE_FES2004_S1 21
#define X_TRACK_REF_TIDE_FES2004_S2 22
#define X_TRACK_REF_TIDE_ANALYSIS 23
#define X_TRACK_REF_TIDE_ANALYSIS_S1 24
#define X_TRACK_REF_TIDE_ANALYSIS_S2 25
#define X_TRACK_REF_NUMVAL 24

#define X_TRACK_SERIE_T_SSHA 0
#define X_TRACK_SERIE_T_DAC_GLOBAL 1
#define X_TRACK_SERIE_T_DAC_GLOBAL_S1 2
#define X_TRACK_SERIE_T_DAC_GLOBAL_S2 3
#define X_TRACK_SERIE_T_DAC_REGIONAL 4
#define X_TRACK_SERIE_T_DAC_REGIONAL_S1 5
#define X_TRACK_SERIE_T_DAC_REGIONAL_S2 6
#define X_TRACK_SERIE_T_DAC_IB 7
#define X_TRACK_SERIE_T_DAC_IB_S1 8
#define X_TRACK_SERIE_T_DAC_IB_S2 9
#define X_TRACK_SERIE_T_DAC_LSA 10
#define X_TRACK_SERIE_T_DAC_SET 11
#define X_TRACK_SERIE_T_TIDE_REGIONAL 12
#define X_TRACK_SERIE_T_TIDE_REGIONAL_S1 13
#define X_TRACK_SERIE_T_TIDE_REGIONAL_S2 14
#define X_TRACK_SERIE_T_TIDE_GOT4_7 15
#define X_TRACK_SERIE_T_TIDE_GOT4_7_S1 16
#define X_TRACK_SERIE_T_TIDE_GOT4_7_S2 17
#define X_TRACK_SERIE_T_TIDE_FES2004 18
#define X_TRACK_SERIE_T_TIDE_FES2004_S1 19
#define X_TRACK_SERIE_T_TIDE_FES2004_S2 20
#define X_TRACK_SERIE_T_TIDE_ANALYSIS 21
#define X_TRACK_SERIE_T_TIDE_ANALYSIS_S1 22
#define X_TRACK_SERIE_T_TIDE_ANALYSIS_S2 23
#define X_TRACK_SERIE_T_RESIDUAL 29

#define X_TRACK_SERIE_T_NUMVAL 24


typedef struct track_data {
  int num_trace;
  int nbpoints;
  int nbpoints_ah;///<number of points on which the harmonic analysis can be computed : some have too little data so impossible
  int nbcycles;
  float *lon;
  float *lat;
  float *lat_ah;///<latitude of the points on which the harmonic analysis is computed
  float *lon_ah;///<idem longitude
  float *mssh;
  int *cycle;
  int *point;
  double **time;
  float **sla;
  float **tide;
  float **mog2d;
  double **amplitude;///<amplitude of the tidal waves at the points at which the harmonic analysis has been computed, dim = nbondes*nbpoints_ah
  double **phase;///<idem phase des ondes
  double **error;///<idem error variance for the inversion

} track_data;


typedef struct track_data_mgr {
  int nbstations;
  float *lat_ah;//latitude des stations
  float *lon_ah;//idem longitude
  double **amplitude;///<amplitude of the tidal waves at the points at which the harmonic analysis has been computed, dim = nbondes*nbstations
  double **phase;///<idem phase des ondes
  double **error;///<idem error variance for the inversion
  char **wave_name; ///<array of the names of the waves to be studied

} track_data_mgr;


extern void handle_error(int status);
extern track_data *read_data_CTOH(const char *dir_input_alti, char *, int *nb_traces_ctoh, int choix_zone);
extern void free_liste_ctoh(track_data *liste, int nb_traces, int nb_ondes);
extern int formatage_CTOH (int nb_traces, track_data *data, const char *dir_ah, char *);
extern int formatage_DUACS (int nb_traces, track_data *data, char *dir_ah);
extern char  **analyse_harmonique (track_data *data,int nb_traces, int choix_zone, int choix_type, char *dir_tide,char *dir_ah, char *dir_ch);

extern int XTRACK_ref_get_nRecords(const string & input);
extern int XTRACK_ref_get_nRecords_ascii(const string & input);
extern int XTRACK_ref_get_nRecords_netcdf(const string & input);

extern std::vector<int> XTRACK_ref_get_indexRecords(const string & input);
extern std::vector<int> XTRACK_ref_get_indexRecords_ascii(const string & input);
extern std::vector<int> XTRACK_ref_get_indexRecords_netcdf(const string & input);

extern int my_nc_get_var1_float(int ncid,const string & varName, const int point, double &value);

extern void XTRACK_ref_get_header(const string & input, int point, double &x, double &y, int &n, int &status);
extern int  XTRACK_ref_get_header_ascii(const string & input, int point, double &x, double &y, int &n);
extern int  XTRACK_ref_get_header_netcdf(const string & input, int point, double &x, double &y, int &n);
extern void XTRACK_ref_get_AllHeaders(const string & input, std::vector<int>, double *x, double *y, int *nrecords, int &status);


extern int XTRACK_get_meteo(serie_t *currentRecord);
extern int XTRACK_get_tide(serie_t *currentRecord);
extern void XTRACK_ref_build_SLA(serie_t *currentRecord, int use_tide, int use_ib);
extern void XTRACK_ref_write_header(FILE *out,int n);

extern void XTRACK_WriteAsciiRecords(FILE *stream, serie_t *record, int *index, size_t size, const char *);

extern void XTRACK_WriteAsciiRecords_obsolete(FILE *stream, serie_t *record, int *index, size_t size);
extern void XTRACK_WriteAsciiRecords_obsolete(FILE *stream, serie_t *record, int *index, size_t size, const char *);

extern int XTRACK_CheckRecords(serie_t & record, int index, const char *rootname);

extern vector<serie_t> ctoh_load_metadata_ASCII (const char *input, int&);
extern vector<serie_t> ctoh_load_metadata_NetCDF(const string & filename);
extern vector<serie_t> ctoh_load_metadata(const string & input, int& status);

extern int ctoh_save_metadata_NetCDF(const string & filename, vector<serie_t> & track);
extern int ctoh_save_metadata_NetCDF(const string & filename,const string & source, vector<serie_t> & track);

extern serie_t XTRACK_ref_get_record_ascii (const char *input, int point, int& status);
extern serie_t XTRACK_ref_get_record_netcdf(const char *input, int point, int& status);
extern serie_t XTRACK_ref_get_record(const string & input, int point, int& status);


#endif
