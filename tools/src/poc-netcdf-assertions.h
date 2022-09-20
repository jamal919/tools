
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief poc-netcdf assertion macros

Included by netcdf-proto.h

*/
/*----------------------------------------------------------------------------*/

#include <iostream>
#include <string> 
#include <netcdf.h>             //for nc_type,NC_...
#include "functions.h"

#if POC_NETCDF_ASSERTIONS_H == 0
#define POC_NETCDF_ASSERTIONS_H 1

extern void nc_check_error(const int stat, const int line, const char *file,int fatal=0);
extern void nc_check_error(const int stat, const int line, const char *file, const string &format, ...);
extern void nc_check_error(const int stat,const char *file,int line,const string & format,...);

/*----------------------------------------------------------------------------*/
/// calls nc_check_error(const int,const int,const char*,const char*,...) without the need to give "__LINE__,__FILE__,"
/** For example : \code NC_CHKERR_LINE_FILE(status,"error with %s",paramName); \endcode */
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define NC_CHKERR_LINE_FILE(s, args... ) do{if(s){nc_check_error(s,__LINE__,__FILE__ , ##args );}}while(0)
#define __NC_CHKERR_LINE_FILE__ NC_CHKERR_LINE_FILE
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// calls nc_check_error(const int,const char*,const int,const char*,...) without the need to give "__FILE__,__LINE__,"
/** For example : \code NC_CHKERR_BASE_LINE(status,"error with %s",paramName); \endcode */
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define NC_CHKERR_BASE_LINE(s, args... ) do{if(s){nc_check_error(s,__FILE__,__LINE__ , ##args );}}while(0)
#define __NC_CHKERR_BASE_LINE__ NC_CHKERR_BASE_LINE
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
/// traps errors
/** For example :
\code if(i<0)NC_TRAP_ERROR(wexit,NC_ENOTVAR,verbose,"%s not found in %s",varName,fileP); \endcode
\code NC_TRAP_ERROR(return,status,verbose,"some_function(\"%s\",\"%s\",) error",fileP,varName); \endcode
In a void function :
\code NC_TRAP_ERROR(return;,status,1,"func(\""+fileP+"\",\""+varName+"\",) returned %d\n",status); \endcode
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define NC_TRAP_ERROR(trap,e,c, args... ) do{\
  if(c>=0)\
    nc_check_error(e,__FILE__,__LINE__ , ##args );\
  if(e) trap(e);\
  }while(0)
  
#define __NC_TRAP_ERROR__ NC_TRAP_ERROR

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

#endif
