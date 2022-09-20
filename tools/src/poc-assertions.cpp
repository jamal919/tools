
/*******************************************************************************

  T-UGO tools, 2006-2014

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief assertion function definition
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h> //for FILE
#include <stdarg.h> //variadic arrays

#include <string.h>
#include <stdlib.h> //for exit,malloc

#include <fstream>

#include "poc-assertions.h"
#include "functions.h"


/*----------------------------------------------------------------------------*/
/// conditionally reports an error and, if so, conditionally exits.
/**
\date reviewed 22 Jul 2011
\author Damien Allain

\param status condition. Does nothing if false
\param *msg message to print
\param line replace by __LINE__ to get line number
\param file replace by __FILE__ to get source file name
\param fatal If ==1, exit.
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void check_error(int status, const char *msg, const int line, const char *file, int fatal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(status != 0) {
    fprintf(stderr, "\n line %d of %s: %s\n", line, file, msg);
    if(fatal == 1) {
      fprintf(stderr,"exiting\n");
      wexit/*avoid perlReplace.pl*/(-1);
      }
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void check_error(const int status, string msg, const int line, const char *file, int fatal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// conditionally reports an error and, if so, conditionally exits.
/**
This is a wrapper for check_error(int,const char*,const int,const char*,int)
*/
/*----------------------------------------------------------------------------*/
{
  check_error(status,msg.c_str(),line,file,fatal);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void check_error(const int status, const int line, const char *file, const char *fmt, ...)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// conditionally reports an error with a custom message.
/**
\date created 1 Aug 2011
\author Damien Allain

\param status condition. Does nothing if false
\param line replace by __LINE__ to get line number
\param file replace by __FILE__ to get source file name
\param "format,..." printf() arguments to message

\code check_error(status,__LINE__,__FILE__,"error with %s",paramName); \endcode
\sa #CHKERR_LINE_FILE
*/
/*----------------------------------------------------------------------------*/
{
  va_list ap;
  if (status != 0) {
    fprintf(stderr, "line %d of %s: ", line, file);
    va_start(ap, fmt);
    vfprintf(stderr, fmt, ap);
    va_end(ap);
    fprintf(stderr, " (%s)\n", strerror(status));
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void wexit(int status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  EXIT POINT : wrapper for exit(), usefull for DEBUGGING

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
/*------------------------------------------------------------------------------
  just put a breakpoint on the line below */

  exit/*avoid perlReplace.pl*/(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void exitIfNull_(void *p,const char *file,int line)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/**

\brief exits with code -1 if a pointer is null

\date 2011-09-15 Damien Allain : creation

Do no call directly. Call #exitIfNull instead.

\param *p If this is null : exits with code -1. Otherwise : do nothing.
\param file replace this with __FILE__
\param line replace this with __LINE__
*/
/*----------------------------------------------------------------------------*/
{
  if(p!=NULL) return;
  fprintf(stderr,"\n%s:%d: null pointer allocation. Exiting with an error code.\n",file,line);
  wexit/*avoid perlReplace.pl*/(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int file_line_start(FILE *f,const char *path,int line)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  r=fprintf(f,"%s:%d:",path,line);
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int file_line_(FILE *f,const char *path,int line, const string &format, ... )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  r=file_line_start(f,path,line);
  
  va_list ap;
  va_start(ap, format);
  r+=vfprintf(f,format.c_str(),ap);
  va_end(ap);
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int base_line_(FILE *f,const char *path,int line, const string &format, ... )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  const char *base;
  base=strrchr0(path,'/');
  r=file_line_start(f,base,line);
  
  va_list ap;
  va_start(ap, format);
  r+=vfprintf(f,format.c_str(),ap);
  va_end(ap);
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void trap_err_exit_(int status,const char *path,int line, const string &format, ... )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  fflush(stdout);
  
  string on_error;
  
  if(status!=0)
    on_error="on error ";
  
  fprintf(stderr,"\n\nExit "+on_error+"in ");
  
  base_line_(stderr,path,line,"\n");
  
  va_list ap;
  va_start(ap, format);
  vfprintf(stderr,format.c_str(),ap);
  va_end(ap);
  
  wexit(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int write_header_if_empty(FILE *f,const string & header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const long
    position=ftell(f);
  
  if(position<0L)
    /* error code */
    return (int)position;
  
  if(position>0L)
    /* not at the start of the file */
    return 0;
  
  int count;
  count=fprintf(f,header+"\n");
  return count;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class A> void assert_array_template(FILE *f,const char *path,int line,const string & name,const char *format,const A & array,int n,int dm=1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(f==0)
    return;
  
  string header,comments,start,separator,end;
  
  /* python : use with import */
  header="from numpy import array";
  comments="# ",
  start="array([",
  separator=",",
  end="])\n";
  
  write_header_if_empty(f,header);
  
  fprintf(f,comments);
  base_line_(f,path,line,"\n");
  
  int i=0,m=0;
  int col=0;
  
  for(;i<n;i++,m+=dm){
    
    if(i==0)
      col+=fprintf(f,name+"="+start);
    else
      col+=fprintf(f,separator);
    
    if(col>=110){
      fprintf(f,"\n");
      col=0;
      }
    
    col+=fprintf(f,format,array[m]);
    }
  
  fprintf(f,end);
  
  fflush(f);
}


const char *doubleFormat="%.15g";
const char *floatFormat="%.7g";


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void assert_array_(FILE *f,const char *path,int line,const string & name,const double *array,int n,int dm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  assert_array_template(f,path,line,name,doubleFormat,array,n,dm);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void assert_array_(FILE *f,const char *path,int line,const string & name,const float *array,int n,int dm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  assert_array_template(f,path,line,name,floatFormat,array,n,dm);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void assert_array_(FILE *f,const char *path,int line,const string & name,const int *array,int n,int dm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  assert_array_template(f,path,line,name,"%d",array,n,dm);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void assert_array_(FILE *f,const char *path,int line,const string & name,const int64_t *array,int n,int dm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  assert_array_template(f,path,line,name,"%lld",array,n,dm);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void assert_array_(FILE *f,const char *path,int line,const string & name,const size_t *array,int n,int dm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  assert_array_template(f,path,line,name,"%lld",array,n,dm);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void assert_array_vector_template(FILE *f,const char *path,int line,const string & name,const char *format,const vector<T> & array,int n=-1,int dm=1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(n<0)
    n=array.size();
  
  assert_array_template(f,path,line,name,format,array,n,dm);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void assert_array_(FILE *f,const char *path,int line,const string &name,const vector<double> & array,int n,int dm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  assert_array_vector_template(f,path,line,name,doubleFormat,array,n,dm);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void assert_array_(FILE *f,const char *path,int line,const string &name,const vector<int> & array,int n,int dm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  assert_array_vector_template(f,path,line,name,"%d",array,n,dm);
}

