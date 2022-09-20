#if POC_TIME_H == 0
#define POC_TIME_H 1

#include <string>
#include <math.h> // for NAN

#include "poc-time-classes.h"

using namespace std;

/*----------------------------------------------------------------------------*/
/// Not A Date
/**
All fields are as impossible as possible : EVEN YEAR 0 DOES NOT EXIST!
\note Negative months and days are tolerated by POSIX date functions

\sa isad() and isnad()
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
date_t const NADate(0,-1,-1,NAN);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/*----------------------------------------------------------------------------*/
/// CNES reference
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
date_t const CNES0(1950,1,1,0.f);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern int julian_day(const date_t &date);
extern int julian_day(int month, int day, int year);

const int CNES0jd __attribute__((unused)) =julian_day(CNES0);

extern int mjd(int month, int day, int year);
extern void calendary(int njd,date_t *actual);

//extern  double cnes_stime(date_t);
//extern  double cnes_dtime(date_t);
extern double time_elapsed(date_t date1, date_t date2);
extern double cnes_time(date_t actual, const char units);

//extern  void getcnesdate(double t,date_t *actual);
//extern  void getcnesdate(double t,date_t *actual,char);
extern  date_t poctime_getdatecnes(double t, char fmt);
extern  date_t poctime_getdatenasa(double t, char fmt);
//extern  void getnasadate(double t,date_t *actual, char fmt);
extern  date_t poctime_getdate(double t,date_t reference,char fmt, double *cnestime);
extern  char *sgetcnesdate(double t);
extern  char *poctime_sdate_cnes(double t, char unit='\0', const char sep=' ');
extern  char *sgetdate(date_t actual, const char sep=' ');

extern int dpm(int month,int year);
extern int day_in_year(const date_t &date);


extern  char *sgetnasadate(double t);
extern  char *sgetnasadatef_01(float *t);
extern  void sgetnasadatef_02(float *t, char *s);

extern int poctime_convert(double *, int, const char, const char);

extern void print_poctime_scan_date_help(int end_date_mode);
extern int poctime_scan_date(const string &str, date_t *start, date_t *end, int verbose=0);
extern date_t poctime_scan_start_date(const char *date_str);
extern date_t poctime_scan_final_date(const char *date_str);

extern bool poctime_date_valid(date_t *date);

extern int poctime_dpm(int month, int year);
extern double ellapsed_time(date_t start, date_t final, const char units);

extern bool poctime_is_after(date_t after, date_t before);
extern bool poctime_is_before(date_t before, date_t after);
extern bool poctime_is_equal(date_t date0, date_t date1);
extern date_t poctime_latest(date_t date0, date_t date1);
extern date_t poctime_soonest(date_t date0, date_t date1);
extern int poctime_count_steps(date_t date_start, date_t date_final, date_t step);
extern date_t poctime_add_step(date_t date, date_t step, date_t stop, date_t *date_end);

#endif
