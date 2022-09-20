
/*******************************************************************************

  T-UGO tools, 2006-2013

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief grib input function declarations
*/
/*----------------------------------------------------------------------------*/

#if POC_GRIB_H == 0
#define POC_GRIB_H 1

#include "../config.h"

#include <string>

#ifdef HAVE_LIBGRIB_API
#include <grib_api.h> /* for grib_get_error_message() */
#endif

using namespace std;

class poc_var_t;
class poc_global_t;
class poc_grid_data_t;
class date_t;

#define POC_GRIB_TIME_NAME "time"
#define POC_GRIB_FILE_PATH_ATT_NAME "filePath"
#define POC_GRIB_FILE_POS_VATT_NAME "fileSeeks"

extern bool is_grib(const string & path);

#ifdef HAVE_LIBGRIB_API
extern int poc_grib_inq(const string &path,poc_global_t *global,int verbose=1);
extern int poc_grib_get_grid_data(const string &path,const poc_var_t &var,poc_grid_data_t *gd,int verbose=0);

extern int poc_grib_get_vara(const string &filename, const poc_var_t &var,int frame, char *z, int verbose=0);
extern int poc_grib_get_vara(const string &filename, const poc_var_t &var,int frame, unsigned char *z, int verbose=0);
extern int poc_grib_get_vara(const string &filename, const poc_var_t &var,int frame, long int *z, int verbose=0);
extern int poc_grib_get_vara(const string &filename, const poc_var_t &var,int frame, double *z, int verbose=0);

extern int poc_grib_get_vara(const string &filename, const poc_var_t &var,int frame, signed char *z, int verbose=0);
extern int poc_grib_get_vara(const string &filename, const poc_var_t &var,int frame, float *z, int verbose=0);
extern int poc_grib_get_vara(const string &filename, const poc_var_t &var,int frame, short *z, int verbose=0);

extern int poc_grib_gettime(const string &path, date_t *origine, double **time, size_t *nframes, poc_global_t *global=0,int verbose=0);
extern int poc_grib_gettime(const string &path, date_t *origine, double **time, int *nframes, poc_global_t *global=0,int verbose=0);
#endif

#endif
