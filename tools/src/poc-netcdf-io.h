
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief poc-netcdf variable input/output declarations
*/
/*----------------------------------------------------------------------------*/

#if POC_NETCDF_IO_HPP == 0
#define POC_NETCDF_IO_HPP 1

#include "poc-netcdf.hpp"


extern bool poc_scale_data(double *data, int length, double scale, double offset, double *mask, int verbose=0);
extern bool poc_scale_data(float *data, int length, double scale, double offset, float *mask, int verbose=0);
extern bool poc_scale_data(float *data, int length, double scale, double offset, float *mask, float specific, int verbose=0);
extern bool poc_scale_data(short int *data, int length, double scale, double offset, short int *mask, short int specific, int verbose=0);
extern bool poc_scale_data(signed char *data, int length, double scale, double offset, signed char *mask, int verbose=0);

extern int poc_scale_data(double *data, int length,const poc_var_t & info, double *mask=0, int verbose=0);


extern int poc_get_vara(int ncid,const poc_var_t &info,int frame, double *z);
extern int poc_get_vara(int ncid,const poc_var_t &info,int frame, float *z);
extern int poc_get_vara(int ncid,const poc_var_t &info,int frame, short int *z);
extern int poc_get_vara(int ncid,const poc_var_t &info,int frame, char *z);
extern int poc_get_vara(int ncid,const poc_var_t &info,int frame, signed char *z);

extern int poc_get_vara(const string &filename, const poc_var_t &info,int frame, double *z, int verbose=0);
extern int poc_get_vara(const string &filename, const poc_var_t &info,int frame, float *z, int verbose=0);
extern int poc_get_vara(const string &filename, const poc_var_t &info,int frame, short int *z, int verbose=0);
extern int poc_get_vara(const string &filename, const poc_var_t &info,int frame, char *z, int verbose=0);
extern int poc_get_vara(const string &filename, const poc_var_t &info,int frame, signed char *z, int verbose=0);

extern int poc_get_vara(int ncid,const poc_var_t &info,double *z,const string & axes="",const indexes_t & indexes= indexes_t() );
extern int poc_get_vara(int ncid,const poc_var_t &info,float *z,const string & axes="",const indexes_t & indexes= indexes_t() );

extern int poc_get_vara(const string &filename,const poc_var_t &info,double *z,const string & axes="",const indexes_t & indexes= indexes_t() , int verbose=0);
extern int poc_get_vara(const string &filename,const poc_var_t &info,float *z,const string & axes="",const indexes_t & indexes= indexes_t() , int verbose=0);
extern int poc_get_vara(const string &filename,const poc_var_t &info,int8_t *z,const string & axes="",const indexes_t & indexes= indexes_t() , int verbose=0);
extern int poc_get_vara(const string &filename,const poc_var_t &info,char *z,const string & axes="",const indexes_t & indexes= indexes_t() , int verbose=0);


extern int poc_get_cvara(int ncid,const string &aname,const string &gname,int frame,complex<double> *z,int verbose=0);
extern int poc_get_cvara(int ncid,const string &aname,const string &gname,int frame,complex<float> *z,int verbose=0);

extern int poc_get_cvara(const string &path,const string &aname,const string &gname,int frame,complex<double> *z,int verbose=0);
extern int poc_get_cvara(const string &path,const string &aname,const string &gname,int frame,complex<float> *z,int verbose=0);
extern int poc_get_cvara(const string &path,const poc_var_t &aname,const poc_var_t &gname,int frame,complex<double> *z,int verbose=0);


extern int poc_get_var(int ncid,const poc_var_t &info,double *z);

extern int poc_get_var(const string & path,const poc_var_t &info,double *z,int verbose=0);


extern double *poc_unscale_data(const double *data, int length, double scale, double offset, double mask=NC_FILL_DOUBLE);
extern float *poc_unscale_data(const float *data, int length, double scale, double offset, float mask=NC_FILL_FLOAT);
extern int *poc_unscale_data(const int *data, int length, double scale, double offset, int mask=NC_FILL_INT);
extern char *poc_unscale_data(const char *data, int length, double scale, double offset, char mask=NC_FILL_CHAR);
extern signed char *poc_unscale_data(const signed char *data, int length, double scale, double offset, signed char mask=NC_FILL_BYTE);


extern int poc_open_and_may_create(const string &path, int *ncid, int verbose=0);


extern int poc_put_vara(int ncid,const poc_var_t &info, int frame,const double *z, int unscale=0);
extern int poc_put_vara(int ncid,const poc_var_t &info, int frame,const float *z, int unscale=0);
extern int poc_put_vara(int ncid,const poc_var_t &info, int frame,const int *z, int unscale=0);
extern int poc_put_vara(int ncid,const poc_var_t &info, int frame,const char *z, int unscale=0);
extern int poc_put_vara(int ncid,const poc_var_t &info, int frame,const signed char *z, int unscale=0);

extern int poc_put_vara(const string &path,const poc_var_t &info, int frame,const double *z, int verbose=0, int unscale=0);
extern int poc_put_vara(const string &path,const poc_var_t &info, int frame,const float *z, int verbose=0, int unscale=0);
extern int poc_put_vara(const string &path,const poc_var_t &info, int frame,const int *z, int verbose=0, int unscale=0);
extern int poc_put_vara(const string &path,const poc_var_t &info, int frame,const char *z, int verbose=0, int unscale=0);
extern int poc_put_vara(const string &path,const poc_var_t &info, int frame,const signed char *z, int verbose=0, int unscale=0);

extern int poc_put_vara(const string &file,const string &var, int frame, const double *z, int verbose=0);
extern int poc_put_vara(const string &file,const string &var, int frame, const float *z, int verbose=0);
extern int poc_put_vara(const string &file,const string &var, int frame, const int *z, int verbose=0);


extern int poc_put_cvara(int ncid,const poc_var_t &ainfo,const poc_var_t &ginfo, size_t frame, const complex<double> *z);
extern int poc_put_cvara(int ncid,const poc_var_t &ainfo,const poc_var_t &ginfo, size_t frame, const complex<float> *z);

extern int poc_put_cvara(const string &path, const poc_var_t &ainfo,const poc_var_t &ginfo, int frame, const complex<double> *z, int verbose=0);
extern int poc_put_cvara(const string &path, const poc_var_t &ainfo,const poc_var_t &ginfo, int frame, const complex<float> *z, int verbose=0);

extern int poc_put_cvara(const string &file,const string &aname,const string &gname, int frame, const complex<double> *z, int verbose=0);
extern int poc_put_cvara(const string &file,const string &aname,const string &gname, int frame, const complex<float> *z, int verbose=0);

extern int poc_put_cvara(const string & path, const poc_var_t & infoTemplate, int frame, const complex<double> *z, int verbose=0);
extern int poc_put_cvara(const string & path, const poc_var_t & infoTemplate, int frame, const complex<float> *z, int verbose=0);


extern int poc_put_var(int ncid,const poc_var_t &info, const double *z, int unscale=0);
extern int poc_put_var(int ncid,const poc_var_t &info, const float *z, int unscale=0);
extern int poc_put_var(int ncid,const poc_var_t &info, const int *z, int unscale=0);

extern int poc_put_var(const string &path, const poc_var_t &info, const double *z, int verbose=0, int unscale=0);
extern int poc_put_var(const string &path, const poc_var_t &info, const float *z, int verbose=0, int unscale=0);
extern int poc_put_var(const string &path, const poc_var_t &info, const int *z, int verbose=0, int unscale=0);

#endif
