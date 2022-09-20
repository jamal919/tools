
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief NetCDF parsing poc-netcdf definitions

Old note : variables routines
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <string.h>
#include <stdarg.h>

#include "tools-structures.h"
#include "poc-netcdf.def"
#include "netcdf-proto.h"

#include "functions.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void nc_check_error(const int stat, const int line, const char *file,int fatal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// conditionally reports a NetCDF error and conditionally exists.
/**
\date reviewed 22 Jul 2011
\author Damien Allain

\param stat NetCDF error status. Does nothing if NC_NOERR
\param line replace by __LINE__ to get line number
\param file replace by __FILE__ to get source file name
\param fatal If 1:exists. Otherwise(default) continues.

\code nc_check_error(status,__LINE__,__FILE__,0); \endcode
*/
/*----------------------------------------------------------------------------*/
{
  if (stat != NC_NOERR) {
    fprintf(stderr, "line %d of %s: %s\n", line, file, nc_strerror(stat));
    if(fatal==1) TRAP_ERR_EXIT(1,"exiting\n");
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void nc_check_error(const int stat, const int line, const char *file, const string &format, ...)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// conditionally reports a NetCDF error with a custom message.
/**
\date reviewed 22 Jul 2011
\author Damien Allain

\param stat NetCDF error status. Does nothing if NC_NOERR
\param line replace by __LINE__ to get line number
\param file replace by __FILE__ to get source file name
\param "format,..." printf() arguments to message

\code nc_check_error(status,__LINE__,__FILE__,"error with "+paramName); \endcode
\sa #NC_CHKERR_LINE_FILE
*/
/*----------------------------------------------------------------------------*/
{
  va_list ap;
  if (stat != NC_NOERR) {
    fprintf(stderr, "line %d of %s: ", line, file);
    va_start(ap, format);
    vfprintf(stderr, format.c_str(), ap);
    va_end(ap);
    fprintf(stderr, " (%d %s)\n",stat, nc_strerror(stat));
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void nc_check_error(const int stat,const char *file,int line,const string & format,...)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// conditionally reports a NetCDF error with a custom message.
/**
\param stat NetCDF error status. Does nothing if NC_NOERR
\param file replace by __FILE__ to get source file name
\param line replace by __LINE__ to get line number
\param "format,..." printf() arguments to message

\code nc_check_error(status,__FILE__,__LINE__,"error with "+paramName); \endcode
\sa #NC_CHKERR_BASE_LINE
*/
/*----------------------------------------------------------------------------*/
{
  va_list ap;
  if (stat != NC_NOERR) {
    fprintf(stderr,"%s:%d:",strrchr0(file,'/'),line);
    va_start(ap, format);
    vfprintf(stderr, format.c_str(), ap);
    va_end(ap);
    fprintf(stderr, " (%d %s)\n",stat, nc_strerror(stat));
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_info (const char* filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, ncid;
  size_t lengthhp;
/*   cdf_fheader header; */
  int ndimsp,nvarsp,ngattsp,unlimdimidp,dim,var,att,nattsp,dimids[50];
  nc_type xtypep;
  char name[NC_MAX_NAME+1];

  char bdum;
  short sdum;
  int idum;
  float fdum;
  double ddum;
  char text[1024];

  #warning cdf_info is still so buggy it has been (temporarily ?) disabled
  __FILE_LINE__(stderr,"cdf_info is still too buggy : skipping.\n");return 0;

  status=nc_open(filename,0,&ncid);
  if(status != NC_NOERR) goto error;

  status=nc_inq(ncid,&ndimsp,&nvarsp,&ngattsp,&unlimdimidp);
  if(status != NC_NOERR) goto error;
  printf("ncid %d ,ndimsp %d,nvarsp %d,ngattsp %d,unlimdimidp %d \n",
          ncid,ndimsp,nvarsp,ngattsp,unlimdimidp);

  for (dim=0;dim<ndimsp;dim++) {
    status=nc_inq_dim(ncid,dim,name,&lengthhp);
    if(status != NC_NOERR) goto error;
    printf("ncid %d ,dimension %d,name %s, length %d \n",ncid,dim,name,lengthhp);
    }

  for (var=0;var<nvarsp;var++) {
    status=nc_inq_var(ncid,var,name,&xtypep,&ndimsp,dimids,&nattsp);
    if(status != NC_NOERR) goto error;
    printf("var %d, name %s, type %d,ndimsp %d, dim %d, natt %d\n",
            var,name,xtypep,ndimsp,dimids[0], nattsp);
    for (att=0;att<nattsp;att++) {
      status=nc_inq_attname(ncid,var,att,name);
      if(status != NC_NOERR) goto error;
      status=nc_inq_att(ncid,var,name,&xtypep,&lengthhp);
      if(status != NC_NOERR) goto error;
/*    printf("var %d ,att %d,name %s, type %d,length %d\n",var,att,name,xtypep,lengthhp);*/
      switch (xtypep) {
        case NC_CHAR:
        ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1024 !
        status=nc_get_att_text(ncid,var,name,text);
        text[lengthhp]='\0';
        if(status != NC_NOERR) goto error;
        printf("att %d, type %d, length %d, name %s, content: %s\n",att,xtypep,lengthhp,name,text);
        break;

        case NC_BYTE:
        ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
        status=nc_get_att_int(ncid,var,name,&idum);
        if(status != NC_NOERR) goto error;
        printf("att %d, type %d, length %d, name %s, content: %d\n",att,xtypep,lengthhp,name,bdum);
        break;

        case NC_SHORT:
        ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
        status=nc_get_att_int(ncid,var,name,&idum);
        if(status != NC_NOERR) goto error;
        printf("att %d, type %d, length %d, name %s, content: %d\n",att,xtypep,lengthhp,name,sdum);
        break;

        case NC_INT:
        ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
        status=nc_get_att_int(ncid,var,name,&idum);
        if(status != NC_NOERR) goto error;
        printf("att %d, type %d, length %d, name %s, content: %d\n",att,xtypep,lengthhp,name,idum);
        break;

        case NC_FLOAT:
        ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
        status=nc_get_att_float(ncid,var,name,&fdum);
        if(status != NC_NOERR) goto error;
        printf("att %d, type %d, length %d, name %s, content: %f\n",att,xtypep,lengthhp,name,fdum);
        break;

        case NC_DOUBLE:
        ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
        status=nc_get_att_double(ncid,var,name,&ddum);
        if(status != NC_NOERR) goto error;
        printf("att %d, length %d, type %d, name %s, content: %lf\n",att,xtypep,lengthhp,name,ddum);
        break;

        }
      }
    }
/*	NC_BYTE =	1,	 signed 1 byte integer */
/*	NC_CHAR =	2,	 ISO/ASCII character */
/*	NC_SHORT =	3,	 signed 2 byte integer */
/*	NC_INT =	4,	 signed 4 byte integer */
/*	NC_FLOAT =	5,	 single precision floating point number */
/*	NC_DOUBLE =	6	 double precision floating point number */

  for (att=0;att<ngattsp;att++) {
    status=nc_inq_attname(ncid,NC_GLOBAL,att,name);
    if(status != NC_NOERR) goto error;
    status=nc_inq_att(ncid,NC_GLOBAL,name,&xtypep,&lengthhp);
    if(status != NC_NOERR) goto error;
/*  printf("global att %d,name %s, type %d,length %d\n", att,name,xtypep,lengthhp);*/
    switch (xtypep) {
      case NC_CHAR:
      ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1024 !
      status=nc_get_att_text(ncid,NC_GLOBAL,name,text);
      if(status != NC_NOERR) goto error;
      printf("global att %d, att %s content: %s\n",att,name,text);
      break;

      case NC_BYTE:
      ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
      status=nc_get_att_int(ncid,NC_GLOBAL,name,&idum);
      if(status != NC_NOERR) goto error;
      printf("global att %d, att %s content: %d\n",att,name,bdum);
      break;

      case NC_SHORT:
      ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
      status=nc_get_att_int(ncid,NC_GLOBAL,name,&idum);
      if(status != NC_NOERR) goto error;
      printf("global att %d, att %s content: %d\n",att,name,sdum);
      break;

      case NC_INT:
      ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
      status=nc_get_att_int(ncid,NC_GLOBAL,name,&idum);
      if(status != NC_NOERR) goto error;
      printf("global att %d, att %s content: %d\n",att,name,idum);
      break;

      case NC_FLOAT:
      ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
      status=nc_get_att_double(ncid,NC_GLOBAL,name,&ddum);
      if(status != NC_NOERR) goto error;
      printf("global att %d, att %s content: %f\n",att,name,ddum);
      break;

      case NC_DOUBLE:
      ///\bug 2011-11-28 Damien Allain : THIS FAILS TO TAKE INTO ACCOUNT THE FACT THAT lengthhp MAY BE >1 !
      status=nc_get_att_double(ncid,NC_GLOBAL,name,&ddum);
      if(status != NC_NOERR) goto error;
      printf("global att %d, att %s content: %lf\n",att,name,ddum);
      break;
      }
    }

  status=nc_close(ncid);
  if(status != NC_NOERR) goto error;

  return(0);

error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void reset_cdfvar(cdfvar_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  info->id=-1;
  info->ndim=0;
  info->type=(nc_type) 0;
  info->natt=0;

  info->att=NULL;
  info->name=NULL;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_diminfo(const int ncid, const int dimid, cdfdim_t *info, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get the cdfdimid_t informations about a dimension from a NetCDF
/**
\param ncid NetCDF file ID
\param dimid
\param *info pointer to the informations
\param verbose
\returns NC_NOERR on success or the NetCDF error status on failure
*/
/*----------------------------------------------------------------------------*/
{
  int unlimdimid,status;
  
  info->destroy();
  info->id=dimid;
  nc_inq_unlimdim(ncid,&unlimdimid);
  info->isunlimited=unlimdimid==dimid;
  
  exitIfNull(info->name=new char[NC_MAX_NAME+1]);
  
  status=nc_inq_dim(ncid,dimid,info->name,&(info->length));
  if(status==NC_NOERR && verbose) {
    cdf_print_diminfo(*info);printf("\n");
    }
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void cdf_print_diminfo(const cdfdim_t info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  printf("%d:%s=%d",info.id,info.name,info.length);
  if(info.isunlimited)printf(" currently");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_attinfo(const int ncid, const int var, const int att, cdfatt_t *info, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  info->id=att;
  exitIfNull(info->name=new char[NC_MAX_NAME+1]);
  status=nc_inq_attname(ncid,var,att,info->name);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_inq_attname error on attribute %d",att);
  status=nc_inq_att(ncid,var,info->name,&info->type,&info->length);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_inq_att error on attribute %s",info->name);
  switch(info->type) {
    case NC_CHAR:/* ISO/ASCII character */
      exitIfNull(info->data=new char[info->length+1]);
      status=nc_get_att_text(ncid,var,info->name,info->data);
      if(status != NC_NOERR) break;
      info->data[info->length]='\0';
      break;
    
    case NC_FLOAT:
      info->data=new char[info->length*sizeof(float)];
      status=nc_get_att_float(ncid,var,info->name,(float *) info->data);
      break;
    
    case NC_DOUBLE:
      info->data=new char[info->length*sizeof(double)];
      status=nc_get_att_double(ncid,var,info->name,(double *) info->data);
      break;
    
    case NC_INT:
      info->data=new char[info->length*sizeof(int)];
      status=nc_get_att_int(ncid,var,info->name,(int *) info->data);
      break;
    
    case NC_SHORT:
      info->data=new char[info->length*sizeof(short)];
      status=nc_get_att_short(ncid,var,info->name,(short *) info->data);
      break;
    
    case NC_BYTE:/* signed 1 byte integer */
      info->data=new char[info->length*sizeof(signed char)];
      status=nc_get_att_schar(ncid,var,info->name,(signed char *) info->data);
      break;
    }
  
  if(status!=NC_NOERR && verbose)NC_CHKERR_LINE_FILE(status,"nc_get_att_... error on attribute %s",info->name);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_varinfo(const int ncid, const int var, cdfvar_t *info, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k;
  int *dimids;

  reset_cdfvar(info);
  if(var<0){
    status=NC_ENOTVAR;
    goto error;
    }
  info->id=var;

/*------------------------------------------------------------------------------
  inquire variables*/

  exitIfNull(info->name=new char[NC_MAX_NAME+1]);

  info->initdim(0);
  status=nc_inq_varndims(ncid,var,&info->ndim);
  if(status != NC_NOERR) goto error;
  info->initdim();
  exitIfNull(dimids=new int[info->ndim]);
  info->initatt(0);
  status=nc_inq_var(ncid,var,info->name,&info->type,NULL,dimids,&info->natt);
  info->initatt();

  for (k=0;k<info->ndim;k++){
    cdf_diminfo(ncid,dimids[k],&info->dim[k],0);
    }
  delete[]dimids;

  for (k=0;k<info->natt;k++) {
    cdf_attinfo(ncid,var,k,&info->att[k],0);
    }

  if(verbose)cdf_print_varinfo(*info);

  status=0;
  return(status);

error:
  nc_check_error(status,__LINE__,__FILE__);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_varinfo(const int ncid, const char *varname, cdfvar_t *info, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**
\date 2011-09-20 Damien ALLAIN : enabled verbose option
*/
{
  int var,status;
  
  reset_cdfvar(info);
  
  status=nc_inq_varid(ncid,varname,&var);
  if(status != NC_NOERR){
    if(verbose)//this function may only be used to check whether a variable exists
      nc_check_error(status,__LINE__,__FILE__,"error on varname %s",varname);
    goto error_close;
    }
  
  info->id=var;
  
  cdf_varinfo(ncid,var,info,verbose);
  if(status != NC_NOERR){nc_check_error(status,__LINE__,__FILE__);goto error_close;}
  
  return(status);

error_close:
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_varinfo(const char *filename, const char *varname, cdfvar_t *info, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**
\date 2011-09-20 Damien ALLAIN : enabled verbose option
*/
{
  int var,ncid,status;
  
  reset_cdfvar(info);
  
  status=nc_open(filename,NC_NOWRITE,&ncid);
  if(status != NC_NOERR) return -1;

  status=nc_inq_varid(ncid,varname,&var);
  if(status != NC_NOERR){
    if(verbose)//this function may only be used to check whether a variable exists
      nc_check_error(status,__LINE__,__FILE__,"error on varname %s",varname);
    goto error_close;
    }
  
  info->id=var;
  
  status=cdf_varinfo(ncid,var,info,verbose);
  if(status != NC_NOERR){nc_check_error(status,__LINE__,__FILE__);goto error_close;}
  
  return nc_close(ncid);

error_close:
  nc_close(ncid);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_varinfo(const char *filename, int var, cdfvar_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;
  
  reset_cdfvar(info);
  
  status=nc_open(filename,NC_NOWRITE,&ncid);
  if(status != NC_NOERR) return -1;

  info->id=var;
  
  cdf_varinfo(ncid,var,info);
  if(status != NC_NOERR){nc_check_error(status,__LINE__,__FILE__);goto error_close;}
  
  return nc_close(ncid);

error_close:
  nc_close(ncid);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_varinfo(const cdfgbl_t &global, const char *varname, cdfvar_t *info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,status=NC_NOERR;
  for(i=0;i<global.nvarsp && strcmp(global.variable[i].name,varname);i++);
  if(i>=global.nvarsp)
    return NC_ENOTVAR;
  *info=global.variable[i];
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void cdf_print_attinfo(const cdfatt_t info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  printf("%d:%s=(%dx",info.id,info.name,info.length);
  switch(info.type){
#define caseType(x) case x:printf(#x);break
  caseType(NC_NAT);
  caseType(NC_BYTE);
  caseType(NC_CHAR);
  caseType(NC_SHORT);
  caseType(NC_INT);
  caseType(NC_FLOAT);
  caseType(NC_DOUBLE);
#ifdef NC_UBYTE
  caseType(NC_UBYTE);
  caseType(NC_USHORT);
  caseType(NC_UINT);
  caseType(NC_INT64);
  caseType(NC_UINT64);
  caseType(NC_STRING);
#else
#warning UPGRADE YOUR NetCDF : NC_UBYTE NC_USHORT NC_UINT NC_INT64 NC_UINT64 NC_STRING are undefined
#endif
#undef caseType
    }
  printf(")");
  switch(info.type){
#define caseType(x,f,t) case x:for(i=0;i<info.length;i++)printf(#f " ",((t*)info.data)[i]);break
//   caseType(NC_BYTE);
//   caseType(NC_SHORT);
  caseType(NC_INT,%d,int);
  caseType(NC_FLOAT,%g,float);
  caseType(NC_DOUBLE,%g,double);
#ifdef NC_UBYTE
//   caseType(NC_UBYTE);
//   caseType(NC_USHORT);
  caseType(NC_UINT,%u,unsigned int);
//   caseType(NC_INT64);
//   caseType(NC_UINT64);
  case NC_CHAR:
  case NC_STRING:
    printf("\"%s\"",info.data);break;
#endif
#undef caseType
    }
  printf(";");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void cdf_print_varinfo(const cdfvar_t &info)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  printf("%d:",info.id);
  switch(info.type){
#define caseType(x) case x:printf(#x);break
  caseType(NC_NAT);
  caseType(NC_BYTE);
  caseType(NC_CHAR);
  caseType(NC_SHORT);
  caseType(NC_INT);
  caseType(NC_FLOAT);
  caseType(NC_DOUBLE);
#ifdef NC_UBYTE
  caseType(NC_UBYTE);
  caseType(NC_USHORT);
  caseType(NC_UINT);
  caseType(NC_INT64);
  caseType(NC_UINT64);
  caseType(NC_STRING);
#endif
#undef caseType
    }
  printf(" %s(",info.name);
  if(info.dim!=NULL)
    for(i=0;i<info.ndim;i++){
      if(i)printf(",");
      cdf_print_diminfo(info.dim[i]);
      }
  else
    printf("%p",info.dim);
  printf("){");
  if(info.att!=NULL)
    for(i=0;i<info.natt;i++){
      cdf_print_attinfo(info.att[i]);
      }
  else
    printf("%p",info.att);
  printf("}\n");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_globalinfo(int ncid,cdfgbl_t *info,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get the cdfgbl_t informations from a NetCDF
/**
\param ncid NetCDF file ID
\param *info pointer to the informations
\param verbose
\returns NC_NOERR on success or the NetCDF error status on failure
*/
/*----------------------------------------------------------------------------*/
  {
  int status;
  int k,dim,var;

/*------------------------------------------------------------------------------
  inquire number of dimensions, variables, global attributes and existence of unlimited dimension*/
  status=nc_inq(ncid,&info->ndimsp,&info->nvarsp,&info->ngattsp,&info->unlimdimid);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"nc_inq error");

  if (verbose) printf("ncid %d ,ndimsp %d,nvarsp %d,ngattsp %d,unlimdimidp %d \n",
                       ncid,info->ndimsp,info->nvarsp,info->ngattsp,info->unlimdimid);

  ///It inquires dimensions with cdf_diminfo()
  exitIfNull(info->dimension = new cdfdim_t[info->ndimsp]);
  for (dim=0;dim<info->ndimsp;dim++){
    status=cdf_diminfo(ncid,dim,&info->dimension[dim],verbose);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"cdf_diminfo error on dim#%d",dim);
    }

  ///It inquires variables with cdf_varinfo(int,int,cdfvar_t*,int)
  exitIfNull(info->variable  = new cdfvar_t[info->nvarsp]);
  for (var=0;var<info->nvarsp;var++) {
    status=cdf_varinfo(ncid, var, &info->variable[var], verbose);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"cdf_varinfo error on var#%d",var);
    }

  ///It inquires attributes with cdf_attinfo(int,int,int,cdfvar_t*,int)
  var=NC_GLOBAL;
  info->initattribute();
  for (k=0;k<info->ngattsp;k++) {
    status=cdf_attinfo(ncid,var,k,&info->attribute[k],verbose);
    if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,verbose,"cdf_attinfo error on att#%d",k);
    if(verbose){cdf_print_attinfo(info->attribute[k]);printf("\n");}
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_globalinfo(const char *filename,cdfgbl_t *info,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Get the cdfgbl_t informations from a NetCDF file
/**
\param *filename file path
\param *info pointer to the informations
\param verbose
\returns NC_NOERR on success or the NetCDF error status on failure
*/
/*----------------------------------------------------------------------------*/
{
  int ncid,status;
  
  status=nc_open(filename,NC_NOWRITE,&ncid);
  if(status != NC_NOERR) goto error;
  
  ///This is basically a wrapper for cdf_globalinfo(int,cdfgbl_t*,int)
  status=cdf_globalinfo(ncid, info, verbose);
  
  nc_close(ncid);
  
  error:
  if(status != NC_NOERR && verbose>=0) NC_CHKERR_BASE_LINE(status,"error on %s",filename);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_identify(const char *filename, const char *varname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int var,ncid,status;

  var=-1;
  status=nc_open(filename,NC_NOWRITE,&ncid);
  if(status != NC_NOERR) goto error;

  status=nc_inq_varid(ncid,varname,&var);
  if(status != NC_NOERR) goto error;

  status=nc_close(ncid);
  if(status != NC_NOERR) goto error;

  return(var);

error:
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_identify(cdfgbl_t g_info, const char *name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int v,id;

  id=-1;
  for(v=0;v<g_info.nvarsp;v++) {
    if(strcmp(g_info.variable[v].name,name)==0) {
      id=v;
      break;
      }
    }

  return(id);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_varid(const char *filename, const char *varname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int var,ncid,status;

  status=nc_open(filename,NC_NOWRITE,&ncid);
  if(status != NC_NOERR) goto error;

  status=nc_inq_varid(ncid,varname,&var);
  if(status!=0) var=-1;

  status=nc_close(ncid);

  return(var);
error:
  return(-1);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_varid(int ncid, const char *varname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int var,status;

  status=nc_inq_varid(ncid,varname,&var);
  if(status!=0) var=-1;

  return(var);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int cdf_identify_dimension(cdfgbl_t g_info, const char *name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int v,id;

  id=-1;
  for(v=0;v<g_info.ndimsp;v++) {
    if(strcmp(g_info.dimension[v].name,name)==0) {
      id=v;
      break;
      }
    }

  return(id);
}

