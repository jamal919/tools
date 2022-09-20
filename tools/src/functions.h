
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief declaration of usefull functions missing from standard libraries
*/
/*----------------------------------------------------------------------------*/

#ifndef _FUNCTIONS_H
#define _FUNCTIONS_H

#include <vector>
#include <string>
#include <map>
#include <string.h> //for memcpy
#include <sys/types.h> //for int64_t
#include <regex.h> /* for regex_t */
#include <errno.h> // for ENOEXEC
#include <complex>

#include <omp.h>
/* OMP_H is not defined by:
- icc that defines __OMP_H
- gcc5 that defines _OMP_H */
#if ( defined(__OMP_H) || ( defined(_OMP_H) ) && !defined(OMP_H) )
#define OMP_H 1
#endif

#include <sys/time.h> //timeval

#ifndef Boolean
#define Boolean char
#endif

using namespace std;


extern int poc_getrf(int n,double *A,int *pivot,int verbose=0);
extern int poc_getrf(int n,complex<double> *A,int *pivot,int verbose=0);

extern int poc_getrs(int n,int nrhs,const double *A,const int *pivot,double *b,int verbose=0);
extern int poc_getrs(int n,int nrhs,const complex<double> *A,const int *pivot,complex<double> *b,int verbose=0);

extern void fflush_std();
extern int asprintf(string & out, const string &format, ... );
extern int sprintf(ostream & out, const string &format, ... );
extern int fprintf(FILE* stream, const string &format, ... );
extern int printf(const string &format, ... );

extern int printc(const string &format, const complex<double> & z);
extern int printc(const string &format, double z);

extern const char *cr,*el,*cuu1;

extern void do_cuu1(int *l);

extern int launchDebugCmd_(const char *file,int line);
#define launchDebugCmd() launchDebugCmd_(__FILE__,__LINE__)

extern string get_extension(const string & path);
extern string replace_extension(const string & path,const string & extension);
extern int get_file_size(const string & path,size_t *size);

extern bool isLetter(char c);
extern const char *strrchr0(const char *s, int c='/');

extern int replace(string *s,const string & match,const string & replacement="",int count=-1);
extern string replace(const string s,const string & match,const string & replacement="",int count=-1);

extern int strncmp(const string &prefix, const string &s);
extern int strrncmp(const string &s1,const string &s2);
extern int strncasecmp(const string & prefix, const string & str);
extern int strrncasecmp(const string & str, const string & suf);
extern void tolower(string *str);
extern char *poc_strdup(const char *source);
extern char *poc_fortran_strdup(const char *source);
extern string poc_fortran_string(const char *source);
extern string canonicalize_file_name(const string & path, bool testEntry=true);
//
extern char *fct_hms(double p);
extern void fct_hms(double p,int *h,int *m,double *s);

extern double scan_angle(const char *s);

extern void fct_xy2polar(float u,float v,float & ro, float & teta);
extern void fct_xy2polar(double u,double v,double & ro, double & teta);

extern void fct_polar2xy(double ro,double teta,double *u,double *v);

extern complex<float> fct_polar2complex(float ro,float teta);

extern int get_date(char *date,const size_t dateL=32uL);
extern string fct_echo(int argc, char *argv[],const char *format="cmdline.%s");

extern char *check_string(const char *word);

extern string dms(double value,bool no0s=false);

extern int count_token(const char *s, const char *delim=" ");
extern int get_token(const char *s, char **params, int nparams);
extern int get_token(const char *s, char ***params);
extern int get_nth_token(const char *s, int parami, char **param);


extern vector<double> get_tokens(string & line, const char* separator);

extern regex_t poc_regcomp_(const char *regex, int cflags,const char *file,int line);
/*----------------------------------------------------------------------------*/
/// wrapper for regcomp()
/**
Will exit with an error code and message if \c regex does not compile.
\return the pattern buffer storage area \c preg
Other parameters like man:/regcomp
*/
/*----------------------------------------------------------------------------*/
#define poc_regcomp(regex,cflags) poc_regcomp_(regex,cflags,__FILE__,__LINE__)

extern vector<string> string_split(const  string  & theString, const  string  & theDelimiter);
extern vector<string> string_split2(const string & line, const string & separator);


extern vector<string> load_filelist(const char *input);
extern int append_filelist(string input, vector<string> & filelist);

extern time_t timeIsOld();

extern double difftime(const struct timeval &b);
extern double difftime(struct timeval *b);

class difftimer_t{
  string baseLineFunc;
  struct timeval t;
  double *durations;
  int id,n,j;
public:
  difftimer_t(const char *file,int line,const char *func,int n0,int id0=0);
#define FILE_LINE_FUNC __FILE__,__LINE__,__func__
  void operator()();
  void destroy();
  ~difftimer_t(){
    destroy();
    };
  };


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

inline int gettimeofday(struct timeval *tv)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
///overload
{
  return gettimeofday(tv,NULL);
}
extern int cmp(const struct timespec &a,const struct timespec &b);


extern int initialize_OPENMP(int gOPENMP_nCPUsMax, FILE *out);
extern int initialize_OPENMP(int gOPENMP_nCPUsMax=-1, int verbose=1);
extern void print_OPENMP_help(const char *prog_name);

typedef struct {
  int64_t on,idle;
  } cpustats_t;
typedef vector<cpustats_t> cpus_stats_t;
extern void get_cpus_stats(cpus_stats_t * stats);
extern void fsleep(double t);
extern int diff_cpus_stats(const cpus_stats_t & stats,vector<double> *ratios,int verbose=0);
extern int diff_cpus_stats(const cpus_stats_t & stats,vector<double> *ratios,FILE *out);
extern int copyTo_tmpdir(string *path,int verbose=0);
extern int copyTo_tmpdir(char **path,int verbose=0);

extern string isFormat_help(const string & indent,const string & inputDefault="",const string & outputDefault="");
extern bool isFormat(char **argv, int *n, string *input, string *output);

extern char *clean_string(const char *word);
extern int decode_format(const string & format, string *input, string *output);

extern void swapValues_(void *a,void *b,size_t n);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> void swapValues(T * a, T * b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* better than a macro to ensure type consistence
will not work on (extremely rare) types declared within functions : use swapValues_() directly instead */
{
  swapValues_(a,b,sizeof(*a));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T>  void swapval(T & a, T & b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  T swap;

  swap = b;
  b = a;
  a = swap;
}

extern bool expire_(const char *file,int line,int ymde, int ymdr=0) __attribute__((warning("the area calling this function has expired or will expire soon")));
/*----------------------------------------------------------------------------*/
/// calls expire_(int,int,int,int) without the need to give "__LINE__,__FILE__,"
/** For example : \code static bool expired=expire(20150420,20160420);
 \endcode */
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define expire(expiry, args... ) expire_(__FILE__,__LINE__,expiry , ##args )
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern int cmp(const struct timeval &a,const struct timeval &b);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <class T> int cmp(T a,T b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(a<b)return -1;
  if(a>b)return 1;
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int cmp(const void *a,const void *b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///ideal for qsort()
/*----------------------------------------------------------------------------*/
{
  return cmp(*(const T*)a,*(const T*)b);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <class T> int is_equal(T a,T b, T precision)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return( fabs(a-b) < precision );
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class T> bool is_between(T a,T b,T c)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  return (a<=b && b<=c)||(a>=b && b>=c);
}


class paire_t {
public:
  int value[2];

  paire_t () {
    value[0]=-1;
    value[1]=-1;
    }
  
  paire_t (int a,int b) {
   value[0]=a;
   value[1]=b;
   }
};


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

struct cmp_str

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/// give char * proper comparisons to map< const char*,> types
{
  bool operator()(char const *a, char const *b){
    return strcmp(a, b) < 0;
    }
};


typedef std::map<const char *, int, cmp_str> toponyms_t;

extern string get_key_list(const toponyms_t & l);
extern string get_key_list(const map<string,int> & l);

#endif
