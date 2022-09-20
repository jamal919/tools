
#if MGR_CONVERTER_H == 0
#define MGR_CONVERTER_H 1

#include "mgr.h"

#define FORMAT_CODE_RMN        1
#define FORMAT_CODE_PAB        2
#define FORMAT_CODE_GNU        3
#define FORMAT_CODE_GLOSS      4
#define FORMAT_CODE_SONEL_HR   5
#define FORMAT_CODE_RADAR      6
#define FORMAT_CODE_HAWAII     7
#define FORMAT_CODE_HAWAII_NC  27
#define FORMAT_CODE_FOREMAN    8
#define FORMAT_CODE_PROFILERS  9
#define FORMAT_CODE_MINH       10
#define FORMAT_CODE_LIST2      11
#define FORMAT_CODE_DMY        12
#define FORMAT_CODE_YMD        13
#define FORMAT_CODE_LEGEND     14
#define FORMAT_CODE_BODC       15
#define FORMAT_CODE_BODC2      16
#define FORMAT_CODE_PUERTOS    17
#define FORMAT_CODE_DART       18
#define FORMAT_CODE_SEINE      19
#define FORMAT_CODE_REFMAR     20
#define FORMAT_CODE_RADAR_RAW  21
#define FORMAT_CODE_NETCDF     22
#define FORMAT_CODE_CANADA     23
#define FORMAT_CODE_LIST       24

#define FORMAT_CODE_USER       99

#define FORMAT_NAME_RMN        "RMN"
#define FORMAT_NAME_PAB        "PAB"
#define FORMAT_NAME_GNU        "GNU"
#define FORMAT_NAME_GLOSS      "GLOSS"
#define FORMAT_NAME_SONEL_HR   "SONEL_HR"
#define FORMAT_NAME_RADAR      "RADAR"
#define FORMAT_NAME_HAWAII     "HAWAII"
#define FORMAT_NAME_HAWAII_NC  "HAWAII_NC"
#define FORMAT_NAME_FOREMAN    "FOREMAN"
#define FORMAT_NAME_PROFILERS  "PROFILERS"
#define FORMAT_NAME_MINH       "MINH"
#define FORMAT_NAME_LIST2      "LIST2"
#define FORMAT_NAME_DMY        "DMY"
#define FORMAT_NAME_YMD        "YMD"
#define FORMAT_NAME_LEGEND     "LEGEND"
#define FORMAT_NAME_BODC       "BODC"
#define FORMAT_NAME_BODC2      "BODC2"
#define FORMAT_NAME_PUERTOS    "PUERTOS"
#define FORMAT_NAME_DART       "DART"
#define FORMAT_NAME_SEINE      "SEINE"
#define FORMAT_NAME_REFMAR     "REFMAR"
#define FORMAT_NAME_RADAR_RAW  "RADAR_RAW"
#define FORMAT_NAME_NETCDF     "NETCDF"
#define FORMAT_NAME_CANADA     "CANADA"
#define FORMAT_NAME_USER       "USER"
#define FORMAT_NAME_LIST       "LIST"

extern  int mgr_init_formats(void);
extern  mooring_t mgr_decode_mooring_info(const string header);
extern int mgr_line_count(string filename);
extern void mgr_list_save_serie(string output, double lon, double lat, double *t, double *sealevel, int nmes);
extern void mgr_list2_save_serie(const string & output,const tseries_t & serie);

extern bool mgr_RMNinfo(string filename, double &lon, double &lat);

extern int mgr_load_serie_nD(const char *filename, int & ncol, double* & time, double** & z, char unit, char *header_format, mooring_t *mooring);

extern int mgr_ParseFormat(const char *filename, const char *format_in);
extern int mgr_load_timeserie(const char *filename, tseries_t **serie, char unit, const char *format_in, int & ncolumns, const char *header_format, mooring_t *mooring, int *status=0);
extern int mgr_decode_position(const char *filename, const char *format_in, const char *header_format, mooring_t *mooring);

extern int mgr_saveforeman(const char *filename,const tseries_t & serie, char unit);
extern int mgr_save_timeserie_gnu(const char *filename,const tseries_t & serie,const mooring_t & mooring, char unit, bool gnu_safe, string time_format="", int verbose=0);
extern int mgr_save_timeserie_NetCDF(const string & filename,const tseries_t *serie, int nseries);
extern int mgr_save_timeserie(const string & rootname,const mooring_t & mooring, const tseries_t & serie, char unit,const string & format, bool gnu_nice, string time_format="", int verbose=0);

extern void mgr_print_formats();
extern int mgr_parse_mooring_header(const char *header_format,const char *line, mooring_t *mooring);
string mgr_build_out_prefix(string f_in, string out_prefix, string *path, string *basename);
extern void mgr_get_path_basename(string in, string *path, string *base);

extern int mgr_converter(const char *file_in, const char *format_in, string delimiter, string time_format, int ncolumns_GNU,
  const char *header_format, const char *mooring_info,
  date_t date_start, date_t date_final, date_t step, int reconstruct,double shift,
  const char *file_out_prefix,const char *format_out, tseries_t ***out_series, double resampling=NAN, bool gnu_nice=1);

extern int mgr_remove_masked_data( tseries_t **series, int nb_series);
extern int mgr_reconstruct(tseries_t *serie,const mooring_t *mooring);

#endif
