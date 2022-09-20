
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\brief assertion macros and related function declarations
*/
/*----------------------------------------------------------------------------*/

#if POC_ASSERTIONS_H == 0
#define POC_ASSERTIONS_H 1

#include <string>
#include <vector>
using namespace std;
#include <stdio.h>

#define FUNCTIONS_H_includes_POC_ASSERTIONS_H 0
#if FUNCTIONS_H_includes_POC_ASSERTIONS_H
extern const char *strrchr0(const char *s, int c);
#else
#include "functions.h" /* for strrchr0() */
#endif

extern void check_error(int status, const char *msg, const int line, const char *file, int fatal=0);
extern void check_error(const int status, string msg, const int line, const char *file, int fatal);
extern void check_error(const int status, const int line, const char *file, const char *fmt, ...);
/*----------------------------------------------------------------------------*/
/// calls check_error(const int,const int,const char*,const char*,...) without the need to give "__LINE__,__FILE__,"
/** For example :
\code __CHKERR_LINE_FILE__(status,"error with %s",paramName); \endcode
\code __CHKERR_LINE_FILE__(ENOEXEC,"not finished");wexit(ENOEXEC); \endcode
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define CHKERR_LINE_FILE(s, args... ) check_error(s,__LINE__,__FILE__ , ##args )
#define __CHKERR_LINE_FILE__ CHKERR_LINE_FILE
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern void wexit(int status);
extern void exitIfNull_(void* p, const char* file, int line);
/*----------------------------------------------------------------------------*/
/** \brief exits with code -1 if a pointer is null

\date 2011-09-15 Damien Allain : creation

This calls nicely exitIfNull_().

\param *p If this is null : exits with code -1. Otherwise : do nothing.

It will print a small message telling which line in the source this occurred.

For example : \code exitIfNull(a=new int[largeInteger]); \endcode
*/
/*----------------------------------------------------------------------------*/
#define exitIfNull(p) exitIfNull_(p,__FILE__,__LINE__)

extern int file_line_(FILE *f,const char *path,int line, const string &format, ... );
extern int base_line_(FILE *f,const char *path,int line, const string &format, ... );
extern void trap_err_exit_(int status,const char *path,int line, const string &format, ... );

/*----------------------------------------------------------------------------*/
/// outputs source file line number to stderr
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define STDERR_LINE(s, args... ) fprintf(stderr,"%d:" s,__LINE__ , ##args )
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// outputs source file line number to stdout
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define STDOUT_LINE(s, args... ) fprintf(stdout,"%d:" s,__LINE__ , ##args )
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// outputs source file BASE name and line number to given file descriptor
/**
\param f file descriptor

Useful for assertions. \sa #__FILE_LINE__
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* To develop, see info:/cpp/Variadic%20Macros */
/* see also info:/cpp-4.7/Swallowing%20the%20Semicolon */
#define BASE_LINE(f,s, args... ) base_line_(f,__FILE__,__LINE__,s , ##args )
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// outputs source file base name and line number to stderr
/**
wrapper for #BASE_LINE */
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define STDERR_BASE_LINE(s, args... ) BASE_LINE(stderr,s , ##args )
/** \todo remove this ubiquitous alias */
#define __ERR_BASE_LINE__ STDERR_BASE_LINE
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// outputs source file base name and line number to stdout
/**
wrapper for #BASE_LINE */
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define STDOUT_BASE_LINE(s, args... ) BASE_LINE(stdout,s , ##args )
/** \todo remove this ubiquitous alias */
#define __OUT_BASE_LINE__ STDOUT_BASE_LINE
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// outputs source file base name, line number and function name to stderr
/**
wrapper for #BASE_LINE */
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define STDERR_BASE_LINE_FUNC(s, args... ) BASE_LINE(stderr," in function %s\n" s,__func__ , ##args )
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// outputs source file base name, line number and function name to stdout
/**
wrapper for #BASE_LINE */
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define STDOUT_BASE_LINE_FUNC(s, args... ) BASE_LINE(stdout," in function %s\n" s,__func__ , ##args )
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// output source file base name and line number to stderr and call return
/**
wrapper for #BASE_LINE
\note This is compatible with void functions
\sa #TRAP_ERR_EXIT
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define TRAP_ERR_RETURN(v,c, args... ) do{if(c>=0){fflush(stdout);BASE_LINE(stderr , ##args );}return v;}while(0)
/** \todo remove this ubiquitous alias */
#define __TRAP_ERR_RETURN__ TRAP_ERR_RETURN
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// output source file base name and line number to stderr and call wexit()
/**
wrapper for #BASE_LINE
\sa #TRAP_ERR_RETURN
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define TRAP_ERR_EXIT(status, args... ) trap_err_exit_(status,__FILE__,__LINE__ , ##args )
/** \todo remove this ubiquitous alias */
#define __TRAP_ERR_EXIT__ TRAP_ERR_EXIT
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// outputs source file COMPILATION PATH name and line number to given file descriptor
/**
\param f file descriptor

Useful for assertions. \sa #BASE_LINE
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define FILE_LINE(f,s, args... ) file_line_(f,__FILE__,__LINE__,s , ##args )
/** \todo remove this ubiquitous alias */
#define __FILE_LINE__ FILE_LINE
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


extern int write_header_if_empty(FILE *f,const string & header);

extern void assert_array_(FILE *f,const char *path,int line,const string &name,const double *array,int n,int dm=1);
extern void assert_array_(FILE *f,const char *path,int line,const string &name,const float  *array,int n,int dm=1);
extern void assert_array_(FILE *f,const char *path,int line,const string &name,const int    *array,int n,int dm=1);
extern void assert_array_(FILE *f,const char *path,int line,const string &name,const int64_t *array,int n,int dm=1);
extern void assert_array_(FILE *f,const char *path,int line,const string &name,const size_t  *array,int n,int dm=1);

extern void assert_array_(FILE *f,const char *path,int line,const string &name,const vector<double> & array,int n=-1,int dm=1);
extern void assert_array_(FILE *f,const char *path,int line,const string &name,const vector<int   > & array,int n=-1,int dm=1);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

#define ASSERT_ARRAY(f,name,array, args... ) assert_array_(f,__FILE__,__LINE__,name,array , ##args )

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


#endif
