
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief definition of usefull functions missing from standard libraries
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h> //for FILE
#include <stdarg.h> //variadic arrays

#include <netcdf.h> //for NC_NOERR
#include <string.h>
#include <stdlib.h> //for exit,malloc

#include <libgen.h> /* for dirname() */

#include <unistd.h> /* for geteuid() */
#include <pwd.h> /* for getpwuid() */

#include <fstream>

#include "poc-assertions.h"
#include "constants.h"

#include "functions.h"

#include "version-macros.def" //for VERSION and REVISION


/* totally fails. More work to be done on this... */
#if 0
bool patchEnv(){
  setenv("OPENBLAS_NUM_THREADS","1",1);
  return true;
}

bool envPatched=patchEnv();
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fflush_std()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  fflush(stderr);
  fflush(stdout);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int asprintf(string & out, const string &format, ... )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  va_list ap;
  char *s;
  
  va_start(ap, format);
  r=vasprintf(&s,format.c_str(),ap);
  va_end(ap);
  
  out+=s;
  free(s);
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sprintf(ostream & out, const string &format, ... )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  va_list ap;
  char *s;
  
  va_start(ap, format);
  r=vasprintf(&s,format.c_str(),ap);
  va_end(ap);
  
  out<<s;
  free(s);
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fprintf(FILE* stream, const string &format, ... )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  va_list ap;
  va_start(ap, format);
  r=vfprintf(stream,format.c_str(),ap);
  va_end(ap);
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int printf(const string &format, ... )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  va_list ap;
  va_start(ap, format);
  r=vprintf(format.c_str(),ap);
  va_end(ap);
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int printc(const string &format, const complex<double> & z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  r=printf("%+"+format+"%+"+format+"j",real(z),imag(z));
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int printc(const string &format, double z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  r=printf("%+"+format,z);
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* see man:/terminfo(5) and output from `infocmp' about the values below */
const char
  *cr   ="\r",       /*< carriage return */      
  *el   ="\x1B[K",   /*< clear to end of line */
  *cuu1 ="\x1B[A";   /*< cursor up one line */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void do_cuu1(int *l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  for(;*l>0;(*l)--)
    fprintf(stderr,"%s",cuu1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void swapValues_(void *a,void *b,size_t n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  void *buffer=malloc(n);
  memcpy(buffer,a,n);
  memcpy(a,b,n);
  memcpy(b,buffer,n);
  free(buffer);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string canonicalize_file_name(const string & path, bool testEntry)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// overloading canonicalize_file_name
/*----------------------------------------------------------------------------*/
{
  char *s;
  
  s=realpath(path.c_str(),0);
  
  if(s==NULL){
    
    if(testEntry)
      return string();
    
    char *ps;
    
    ps=poc_strdup(path.c_str());
    
    string r;
    r=dirname(ps);
    delete[]ps;
    
    r=canonicalize_file_name(r);
    r+="/";
    
    ps=poc_strdup(path.c_str());
    r+=basename(ps);
    delete[]ps;
    
    return r;
    }
  
  string r(s);
  free(s);
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int get_date(char *date,const size_t dateL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Puts the current date and time in a string
/**
\param date string
\param dateL maximum string length
\returns 0 on success or errno if error
*/
/*----------------------------------------------------------------------------*/
{
  time_t call_time;
  struct tm call_tm;
  call_time=time(NULL);
  size_t fitted;
  if(call_time==-1){
    CHKERR_LINE_FILE(errno,"time() error");
    return errno;
    }
  if(localtime_r(&call_time,&call_tm)==NULL){
    CHKERR_LINE_FILE(errno,"localtime_r() error");
    return errno;
    }
  /**This puts the date and time as yyyy-mm-dd HH:MM:SS TZ format, as shown below :
  \code /**/ // COMPILED CODE BELOW !!!
  fitted=strftime(date,dateL,"%F %T %Z",&call_tm); /** \endcode */ //returns the number of characters placed in the array s, not including the terminating null byte, provided the string, including the terminating null byte, fits. Otherwise, it returns 0, and the contents of the array is undefined.
  if(fitted==0){
    CHKERR_LINE_FILE(ERANGE,"strftime() error");
    return ERANGE;/* Result too large */
    }
  return 0;//errno is never set to zero by any system call or library function
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string fct_echo(int argc, char *argv[],const char *format)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Records the date and time of a program launch, as well as its arguments
/** in a file on the current directory.

The arguments are the same as those of the main function
\param argc
\param argv
\returns the command line
For example: \code int main( int argc, char*[] argv ){
  ...
  fct_echo( argc, argv);
  ...
  }
  \endcode
*/
/*----------------------------------------------------------------------------*/
{
  string s;
  int n,pos,status;
  FILE *echo=NULL;
  const int
    dateL=32, /* should only need 32 chars for storing date + time + timezone */
    hostL=65; /* as per man:/gethostname(2), see definition of utsname.nodename in <sys/utsname.h> */
  char date[dateL],host[hostL],*path;
  struct passwd *pw;
  const uid_t uid=geteuid();

  ///-get the date and time with get_date()
  n=get_date(date,dateL);
  if(n!=0){
    CHKERR_LINE_FILE(n,"get_date() error");
    return NULL;
    }

  s=argv[0];
  pos=s.rfind("/")+1;
  asprintf(&path,format,s.substr(pos).c_str());

  for(n=1;n < argc;n++) {
    char *argvn=argv[n];
    string q;
    
    if(strpbrk(argvn," '<>#;()*?")!=0 or argvn[0]=='\0'){
      q="\"";
      }
    else if( strpbrk(argvn,"$\"")!=0 ){
      q="'";
      }
    
    s+=" "+q+argvn+q;
    }

  echo=fopen(path,"a");
  if(echo==0) {
    CHKERR_LINE_FILE(errno,"fopen error on %s",path);
    goto error;
    }
  fprintf(echo,"\n%s",date);
  pw=getpwuid(uid);
  status=gethostname(host,hostL);
  fprintf(echo," %s@%s",pw->pw_name,host); /* user@machine */
  fprintf(echo," " HG_REV); /* Mercurial revison */
  fprintf(echo,"\n%s\n",s.c_str());
  fclose(echo);

error:
  free(path);
  return s;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fct_hms(double p,int *h,int *m,double *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(isnan(p)){
    *m=-1;
    *s=NAN;
    return;
    }
  
  double t;

  *h=floor(p);
  t=(p-*h)*60.0;
  *m=floor(t);
  *s=(t-*m)*60.0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *fct_hms(double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int h,m;
  char *s=NULL;
  double t;

  fct_hms(p,&h,&m,&t);
  
  s=new char[10];
  
  sprintf(s,"%2.2dh %2.2dm %6.3lfs",h,m,t);
  
  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double fct_fromhms(int d, int m, double s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double t;

  t=(double) d+(double) m/60.+s/3600.;

  return(t);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double scan_angle(const char *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double r;
  char *unit;

  r=strtod(s,&unit);
  if(unit==&s[0])
    return NAN;
  switch(tolower(*unit)){
  case '"':
  case 's':
    r/=60.;
  case '\'':
  case 'm':
    r/=60.;
  case '\0':
  case 'd':
    break;
  case 'r':
    r*=r2d;
    break;
  case 'g':
    r*=200./M_PI;
    break;
  default:
    STDERR_BASE_LINE_FUNC("unit %s not recognised\n",unit);
    return NAN;
    }

  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> void fct_xy2polar_template(T u, T v, T & ro, T & teta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  teta en degre */
  if((u != 0.0) || (v != 0.0)) {
    teta=atan2(v,u);
    teta=teta *r2d;
    if(teta < 0.0) teta+=360.;
    }
  else
    teta=0.0;

  ro=sqrt(u*u + v*v);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fct_xy2polar(double u, double v, double & ro, double & teta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  fct_xy2polar_template(u, v, ro, teta);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fct_xy2polar(float u, float v, float & ro, float & teta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  fct_xy2polar_template(u, v, ro, teta);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fct_polar2xy(double ro,double teta,double & u,double & v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------
  teta en degre */
  u=ro*cos(teta*d2r);
  v=ro*sin(teta*d2r);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<float> fct_polar2complex(float ro, float teta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<float> result;
  
  result=polar<float>(ro,teta*d2r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  complex<double> fct_polar2complex(double ro, double teta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> result;
  
  result=polar<double>(ro,teta*d2r);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string dms(double value,bool no0s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///converts a number to degree, minute and second string
/**
\param value
\return the string
*/
/*----------------------------------------------------------------------------*/
{
  string str;
  double d,m,factor=1.;

/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  if(value<0){
    value*=-1.;
    factor=-1.;
    }
  
  d=floor(value);
  value-=d;
  if(1.-value<1./3600e3){
    d++;
    value=0.;
    }
  
  if(d!=0) {
    asprintf(str,no0s?"%d°":"%4d°",(int) (d*factor));
    }
  else if(factor==1) {
    asprintf(str,no0s?"%d°":"%4d°",(int) d);
    }
  else
    {
    asprintf(str,no0s?"%d°":"  -%1d°",(int) d);
    }

  if(no0s && (value==0.0) ) goto end;
  
  value*=60.0;
  m=floor(value);
  value-=m;
  if(1.-value<1./3600e3){
    m++;
    value=0.;
    }
  
  asprintf(str,"%02d'",(int)m);
  value*=60.0;
  if(no0s && abs(value)<1e-6) goto end;
  
  {
  int i;
  string word;
  /* NOTE: 1" is ~31m */
  asprintf(word,"%05.2f",value);
  if(no0s){
    for(i=word.size()-1;i>1 && strchr(".0",word[i]);i--);
    word[i+1]='\0';
    }
  str+=word+'"';
  }
  
  end:
  return str; /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<string> load_filelist(const char *input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<string> filelist;
  ifstream in(input);
  string filename;

  if(!in){
    STDERR_BASE_LINE("unable to open \"%s\" : ",input);perror("");
    return(filelist);
    }

  while(in.good()){
    getline(in,filename);
    if(in.eof()){
      break;
      }
    filelist.push_back(filename);
    };

  in.close();
  
  return(filelist);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  vector<string> load_filelist(string input)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<string> filelist;
  
  filelist=load_filelist(input.c_str());
  
  return(filelist);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int append_filelist(string input, vector<string> & filelist)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ifstream in(input.c_str());
  string filename;

  if(!in){
    STDERR_BASE_LINE("unable to open \"%s\" : ",input.c_str());perror("");
    return(-1);
    }

  while(in.good()){
    getline(in,filename);
    if(in.eof()){
      break;
      }
    filelist.push_back(filename);
    };

  in.close();
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void deletep_void_void(void **p,void (*F)(void *))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// delete a void pointer and set it to \c NULL
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  if(!*p) return;
  
  F(*p);
  
  *p=NULL;
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> string get_key_list_template(const T & l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string keys;
  
  for(typename T::const_iterator it=l.begin(); it!=l.end(); it++) {
    if(keys!="")
      keys+=" ";
    keys+=it->first;
    }
  
  return keys;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string get_key_list(const toponyms_t & l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string keys;
  
  keys=get_key_list_template(l);
  
  return keys;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  string get_key_list(const map<string,int> & l)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string keys;
  
  keys=get_key_list_template(l);
  
  return keys;
}
