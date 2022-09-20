
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief NetCDF variables definition routines

BECAUSE OF SOME FUNCTIONS, IT NEEDS TO BE COMBINED WITH geo-02.cpp (libgeo.a)
*/
/*----------------------------------------------------------------------------*/

#include "version-macros.def"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <netcdf.h>

#include "tools-structures.h"
#include "poc-netcdf.def"
#include "netcdf-proto.h"
#include "geo.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_createfile(const char *filename, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
   int  ncid; /* netCDF id */
   char *value=NULL;
   int status, lgth;
//   int cmode=NC_CLOBBER | NC_NETCDF4|NC_CLASSIC_MODEL;
   int cmode;
   
   if(mode==-1) cmode=NC_CLOBBER;
   else cmode=mode;

   /* enter define mode */
   status = nc_create(filename, cmode, &ncid);
   nc_check_error(status,__LINE__,__FILE__);

   /* define dimensions */

   /* define variables */

   /* assign global attributes */
   if(value!=NULL) free(value);

   value=strdup("CF1.0, POC revision 2010");
   lgth=strlen(value);
   status = nc_put_att_text(ncid, NC_GLOBAL, "Conventions", lgth,value);
   nc_check_error(status,__LINE__,__FILE__);
   if(value!=NULL) free(value);
   value=NULL;

   value=strdup("regular");
   lgth=strlen(value);
   status = nc_put_att_text(ncid, NC_GLOBAL, "grid_type", lgth, "regular");
   nc_check_error(status,__LINE__,__FILE__);
    if(value!=NULL) free(value);
   value=NULL;

   lgth=strlen(filename);
   status = nc_put_att_text(ncid, NC_GLOBAL, "file_name", lgth, filename);
   nc_check_error(status,__LINE__,__FILE__);

   asprintf(&value,"file produced around line %d of " __FILE__ " of " PACKAGE_STRING " " REVISION,__LINE__);
   lgth=strlen(value);
   status = nc_put_att_text(ncid, NC_GLOBAL, "production", lgth, value);
   nc_check_error(status,__LINE__,__FILE__);
   if(value!=NULL) free(value);
   value=NULL;

  time_t creation_time;
  status = time(&creation_time);

  value=strdup(ctime(&creation_time));
  lgth=strlen(value);
  value[lgth-1]=0;
  lgth--;
  status = nc_put_att_text(ncid, NC_GLOBAL, "creation_date", lgth, value);
  nc_check_error(status,__LINE__,__FILE__);
  if(value!=NULL) free(value);
  value=NULL;

   /* leave define mode */
   status = nc_enddef (ncid);
   nc_check_error(status,__LINE__,__FILE__);
   status = nc_close(ncid);
   nc_check_error(status,__LINE__,__FILE__);

   return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_createfile(const char *filename, const cdfgbl_t &global)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,dim_id,ncid;
  char *value=NULL;
  int status, lgth;
  size_t fileSize;

  status=get_file_size(filename, &fileSize);

//   status = nc_open(filename, NC_WRITE, &ncid);
//   nc_check_error(status,__LINE__,__FILE__);
  if(status==-1) {
    status = nc_create(filename, NC_CLOBBER, &ncid);
    nc_check_error(status,__LINE__,__FILE__);
    if(status!=NC_NOERR)return status;
    }
  else {
    status = nc_open(filename, NC_WRITE, &ncid);
    nc_check_error(status,__LINE__,__FILE__);
    if(status!=NC_NOERR)return status;
    status = nc_redef(ncid);
    nc_check_error(status,__LINE__,__FILE__);
    if(status!=NC_NOERR)goto error;
    }

/*------------------------------------------------------------------------------
  create dimensions */
  for(k=0;k<global.ndimsp;k++) {
    status=nc_inq_dimid(ncid,global.dimension[k].name,&dim_id);
    if (status ==  NC_EBADDIM) {
      if(k==global.unlimdimid) {
        status=nc_def_dim(ncid,global.dimension[k].name, NC_UNLIMITED, &dim_id);
        }
      else {
        status=nc_def_dim(ncid,global.dimension[k].name, global.dimension[k].length, &dim_id);
        }
      nc_check_error(status,__LINE__,__FILE__);
      }
    }

/*------------------------------------------------------------------------------
  create variables */
  for(k=0;k<global.nvarsp;k++) {
    status=create_ncvariable(ncid,&global.variable[k]);
    NC_CHKERR_LINE_FILE(status,"create_ncvariable error with %s on %s",global.variable[k].name,filename);
    }


//   var=NC_GLOBAL;
//   info->initattribute();
//   for (k=0;k<info->ngattsp;k++) {
//     info->attribute[k].name=(char *) malloc(1024);
//     statusus=nc_inq_attname(ncid,var,k,info->attribute[k].name);
//     if(statusus != NC_NOERR) goto error;
//     statusus=nc_inq_att(ncid,var,info->attribute[k].name,&info->attribute[k].type,&info->attribute[k].length);
//     if(statusus != NC_NOERR) goto error;
//     switch(info->attribute[k].type) {
//       case NC_CHAR:
//         info->attribute[k].data=(char *) malloc(info->attribute[k].length+1);
//         statusus=nc_get_att_text(ncid,var,info->attribute[k].name,info->attribute[k].data);
//         if(statusus != NC_NOERR) goto error;
//         info->attribute[k].data[info->attribute[k].length]='\0';
//         break;
//
//       case NC_FLOAT:
//         info->attribute[k].data=new char[sizeof(float)];
//         statusus=nc_get_att_float(ncid,var,info->attribute[k].name,(float *) info->attribute[k].data);
//         if(statusus != NC_NOERR) goto error;
//         break;
//
//       case NC_DOUBLE:
//         info->attribute[k].data=new char[sizeof(double)];
//         statusus=nc_get_att_double(ncid,var,info->attribute[k].name,(double *) info->attribute[k].data);
//         if(statusus != NC_NOERR) goto error;
//         break;
//
//       case NC_INT:
//         info->attribute[k].data=new char[sizeof(int)];
//         statusus=nc_get_att_int(ncid,var,info->attribute[k].name,(int *) info->attribute[k].data);
//         if(statusus != NC_NOERR) goto error;
//         break;
//       }
//     }

  asprintf(&value,"file produced around line %d of " __FILE__ " of " PACKAGE_STRING " " REVISION,__LINE__);
  lgth=strlen(value);
  status = nc_put_att_text(ncid, NC_GLOBAL, "production", lgth, value);
  nc_check_error(status,__LINE__,__FILE__);
  if(value!=NULL) free(value);
  value=NULL;

  time_t creation_time;
  status = time(&creation_time);

  value=strdup(ctime(&creation_time));
  lgth=strlen(value);
  value[lgth-1]=0;
  lgth--;
  status = nc_put_att_text(ncid, NC_GLOBAL, "creation_date", lgth, value);
  nc_check_error(status,__LINE__,__FILE__);
  if(value!=NULL) free(value);
  value=NULL;


  /* leave define mode */
  status = nc_enddef (ncid);
  NC_CHKERR_LINE_FILE(status,"nc_enddef error on %s",filename);
error:
  nc_close(ncid);

  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

char *poc_getdate(date_t actual)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double ftime,fhours,fminutes,fseconds;
  int    hours,minutes,seconds;
  char *s=NULL;

  ftime=(double) actual.second;

  fhours=floor(ftime/3600.);

  fminutes=floor((ftime-fhours*3600.)/60.);

  fseconds=ftime-fhours*3600.-fminutes*60.;

  hours=int(floor(fhours+0.5));
  minutes=int(floor(fminutes+0.5));
  seconds=int(floor(fseconds+0.5));

  exitIfNull(
    s=(char *)malloc(256*sizeof(char))
    );
  sprintf(s,"%4.4d-%s-%2.2d %2.2d:%2.2d:%2.2d",actual.year,uc_names[actual.month-1],actual.day,hours,minutes,seconds);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int create_ncvariable(int ncid, cdfvar_t *variable, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  k,len;

  int status;

  int *dimids;

  cdfatt_t *att;

  dimids=new int[variable->ndim];
  for (k=0;k<variable->ndim;k++) {
    status=nc_inq_dimid(ncid,variable->dim[k].name,&variable->dim[k].id);
    if(status!=0) {
      printf("nc_inq_dimid status %d: dimension %s length=%d\n",status,variable->dim[k].name,variable->dim[k].length);
      nc_check_error(status,__LINE__,__FILE__);
      }
    dimids[k]=variable->dim[k].id;
    }

/*------------------------------------------------------------------------------
  enter define mode */
  status = nc_def_var(ncid, variable->name, variable->type, variable->ndim, dimids, &variable->id);
  delete[]dimids;
  if(status!=NC_NOERR) {
/*------------------------------------------------------------------------------
    variable already exist, hope it is allright */
    if(status==-42) {
      if(verbose) printf("variable %s already exists, will assume it is ok\n",variable->name);
      return(0);
      }
    NC_TRAP_ERROR(return,status,1,"nc_def_var error with %s",variable->name);
    }

  for (k=0;k<variable->natt;k++){
    att=&variable->att[k];
    switch (att->type) {
      case NC_CHAR:
        len= strlen(att->data);
        status = nc_put_att_text(ncid, variable->id, att->name,len,att->data);
        nc_check_error(status,__LINE__,__FILE__);
        break;

      case NC_BYTE:
        status =-1;
        nc_check_error(status,__LINE__,__FILE__);
        break;

      case NC_SHORT:
        status = nc_put_att_short(ncid,variable->id, att->name, NC_SHORT, att->length, (short *)att->data);
        nc_check_error(status,__LINE__,__FILE__);
        break;

      case NC_INT:
        status = nc_put_att_int(ncid,variable->id, att->name, NC_INT, att->length, (int *)att->data);
        nc_check_error(status,__LINE__,__FILE__);
        break;

      case NC_FLOAT:
        status = nc_put_att_float(ncid,variable->id, att->name, NC_FLOAT, att->length, (float *)att->data);
        nc_check_error(status,__LINE__,__FILE__);
        break;

      case NC_DOUBLE:
        status = nc_put_att_double(ncid,variable->id, att->name, NC_DOUBLE, att->length, (double *)att->data);
        nc_check_error(status,__LINE__,__FILE__);
        break;
      }
    }

   return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cdf_createvariable(const char *filename,cdfvar_t *variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  ncid;	/* netCDF id */
  int  k,len;

  int stat;

  int *dimids;

  /* attribute vectors */
  short short_value;
  float float_value;
  double double_value;

/*------------------------------------------------------------------------------
  enter define mode */
  stat = nc_open(filename, NC_WRITE, &ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat=nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  dimids=new int[variable->ndim];
  for (k=0;k<variable->ndim;k++) {
    stat=nc_inq_dimid(ncid,variable->dim[k].name,&variable->dim[k].id);
    nc_check_error(stat,__LINE__,__FILE__);
    dimids[k]=variable->dim[k].id;
    }

  stat = nc_def_var(ncid, variable->name, variable->type, variable->ndim, dimids, &variable->id);
  delete[]dimids;
  NC_TRAP_ERROR(wexit,stat,1,"nc_def_var((\"%s\"),(\"%s\"),...) error");

  for (k=0;k<variable->natt;k++) {
      switch (variable->att[k].type) {
        case NC_CHAR:
          len= strlen(variable->att[k].data);
          stat = nc_put_att_text(ncid, variable->id, variable->att[k].name,len,variable->att[k].data);
          nc_check_error(stat,__LINE__,__FILE__);
          break;

        case NC_BYTE:
          stat =-1;
          nc_check_error(stat,__LINE__,__FILE__);
          break;

        case NC_SHORT:
          short_value=*((short *) variable->att[k].data);
          stat = nc_put_att_short(ncid, variable->id, variable->att[k].name,NC_SHORT,1,&short_value);
          nc_check_error(stat,__LINE__,__FILE__);
          break;

        case NC_INT:
          stat =-1;
          nc_check_error(stat,__LINE__,__FILE__);
          break;

        case NC_FLOAT:
          float_value=*((float *) variable->att[k].data);
          stat = nc_put_att_float(ncid,variable->id, variable->att[k].name, NC_FLOAT, 1, &float_value);
          nc_check_error(stat,__LINE__,__FILE__);
          break;

        case NC_DOUBLE:
          double_value=*((double *) variable->att[k].data);
          stat = nc_put_att_double(ncid,variable->id, variable->att[k].name, NC_DOUBLE, 1, &double_value);
          nc_check_error(stat,__LINE__,__FILE__);
          break;
        }
    }

  stat = nc_enddef (ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int create_ncvariable(const char *filename,cdfvar_t *variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int ncid,status;

/*------------------------------------------------------------------------------
  open file in write mode */
  status = nc_open(filename, NC_WRITE, &ncid);
  if(status != NC_NOERR) goto error;
  
  status = nc_inq_varid(ncid, variable->name, &(variable->id));
 
  if(status==NC_NOERR) {
    printf("variable %s already exists, will assume it is ok\n",variable->name);
    nc_close(ncid);
    return(0);
    }
  
/*------------------------------------------------------------------------------
  enter define mode */
  status=nc_redef(ncid);
  if(status != NC_NOERR) goto error;

  create_ncvariable(ncid,variable);

  status = nc_enddef (ncid);
  error:
  nc_check_error(status,__LINE__,__FILE__);
  nc_close(ncid);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void add_att(cdfvar_t *v,cdfatt_t a)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
///add attribute to a variable
/**
\param a attribute to append
Note: an attribute of the same name overwrites the previous one
*/
{
  //NOTE: avoid putting breakpoints here !!!
  //use assertions instead
  //__FILE_LINE__("\n");
  cdfatt_t *att0;
  att0=v->att;
  if(v->att==NULL)
    v->natt=0;//fail safe
  v->att=new cdfatt_t[v->natt+1];
  if(v->att!=NULL){
    memcpy(v->att,att0,v->natt*sizeof(cdfatt_t));
    delete[]att0;
    }
  v->att[v->natt]=a;
  //memcpy(&v->att[v->natt],&att0,sizeof(cdfatt_t));
  v->natt++;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int create_ncvariable_indfile(char *filename,cdfvar_t *variable,char *filegrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int  ncid,ncidgrid; /* netCDF id */
  int  k,len;

   int stat;

   int *dimids;

   /* attribute vectors */
   float float_value;
   double double_value;

/*------------------------------------------------------------------------------
   enter define mode */
   stat = nc_open(filename, NC_WRITE, &ncid);
   stat = nc_open(filegrid, NC_NOWRITE, &ncidgrid);
   nc_check_error(stat,__LINE__,__FILE__);
   stat=nc_redef(ncid);
   nc_check_error(stat,__LINE__,__FILE__);

   dimids=new int[variable->ndim];
   for (k=0;k<variable->ndim;k++) {
       stat=nc_inq_dimid(ncid,variable->dim[k].name,&variable->dim[k].id);
       if (stat == NC_EBADDIM) {
         stat=nc_inq_dimid(ncidgrid,variable->dim[k].name,&variable->dim[k].id);
         nc_check_error(stat,__LINE__,__FILE__);
         stat = nc_inq_dimlen(ncidgrid, variable->dim[k].id, &variable->dim[k].length);
         nc_check_error(stat,__LINE__,__FILE__);
         stat = nc_def_dim (ncid,variable->dim[k].name,variable->dim[k].length, &variable->dim[k].id);
         nc_check_error(stat,__LINE__,__FILE__);
       }
       dimids[k]=variable->dim[k].id;
     }
     
   stat = nc_def_var(ncid, variable->name, variable->type, variable->ndim, dimids, &variable->id);
   delete[]dimids;
   nc_check_error(stat,__LINE__,__FILE__);

   for (k=0;k<variable->natt;k++) {
       switch (variable->att[k].type)
         {
         case NC_CHAR:
           len= strlen(variable->att[k].data);
           stat = nc_put_att_text(ncid, variable->id, variable->att[k].name,len,variable->att[k].data);
           nc_check_error(stat,__LINE__,__FILE__);
           break;

         case NC_BYTE:
           stat =-1;
           nc_check_error(stat,__LINE__,__FILE__);
           break;

         case NC_SHORT:
           stat =-1;
           nc_check_error(stat,__LINE__,__FILE__);
           break;

         case NC_INT:
           stat =-1;
           nc_check_error(stat,__LINE__,__FILE__);
           break;

         case NC_FLOAT:
           float_value=*((float *) variable->att[k].data);
           stat = nc_put_att_float(ncid,variable->id, variable->att[k].name, NC_FLOAT, 1, &float_value);
           nc_check_error(stat,__LINE__,__FILE__);
           break;

         case NC_DOUBLE:
           double_value=*((double *) variable->att[k].data);
           stat = nc_put_att_double(ncid,variable->id, variable->att[k].name, NC_DOUBLE, 1, &double_value);
           nc_check_error(stat,__LINE__,__FILE__);
           break;
         }
     }

   stat = nc_enddef (ncid);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_close(ncid);
   nc_check_error(stat,__LINE__,__FILE__);
   stat = nc_close(ncidgrid);
   nc_check_error(stat,__LINE__,__FILE__);

   return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_sphericalgrid_xy(const char *input, const char *name,const grid_t & grid, pocgrd_t *cdfgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  int k,l,status;
  const size_t nd=2;
  char     *dname[nd];
  char     *tmp[nd];
  int dimlgth[nd];
  int ncid,stat,dim_id;

  cdfgrid->id=-1;
  char *lonname=NULL, *latname=NULL, *levelsname=NULL;
  char *xname=NULL, *yname=NULL;

/*------------------------------------------------------------------------------

  Defines a grid to which variable will be fully or partially connected.

  1) defines dimensions
     beware the unlimited dimension

  2) creates standard horizontal and vertical axis


  3) returns informations to be passed to variable definition and creation routines

*------------------------------------------------------------------------*/

  l=strlen(name);
  if(l==0) {
    dname[0]=poc_strdup("nx");
    dname[1]=poc_strdup("ny");
    lonname=poc_strdup("lon");
    latname=poc_strdup("lat");
    xname=poc_strdup("x");
    yname=poc_strdup("y");
    }
  else {
    for(k=0;k<4;k++) {
      exitIfNull(
        dname[k]=new char[l+3]
        );
      }
    sprintf(dname[0],"nx_%s",name);
    sprintf(dname[1],"ny_%s",name);
    exitIfNull(
      lonname=new char[l+5]
      );
    exitIfNull(
      latname=new char[l+5]
      );
    exitIfNull(
      levelsname=new char[l+7]
      );
    sprintf(lonname,"lon_%s",name);
    sprintf(latname,"lat_%s",name);
    }

/*------------------------------------------------------------------------------
  enter define mode */
  stat=nc_open(input, NC_WRITE, &ncid);
  if(stat==ENOENT){
    stat=poc_createfile(input, NC_64BIT_OFFSET);
    nc_check_error(stat,__LINE__,__FILE__);
    stat=nc_open(input, NC_WRITE, &ncid);
    }
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create dimensions */
  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  for(k=0;k<2;k++) {
    stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
    nc_check_error(stat,__LINE__,__FILE__);
    }
  stat = nc_enddef (ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create axis variable */

  exitIfNull(cdfgrid->x   = new cdfvar_t[1]);
  exitIfNull(cdfgrid->y   = new cdfvar_t[1]);
  exitIfNull(cdfgrid->lon = new cdfvar_t[1]);
  exitIfNull(cdfgrid->lat = new cdfvar_t[1]);
  exitIfNull(cdfgrid->z   = new cdfvar_t[1]);

  switch (grid.modeH) {
    case 0:
      
    case 1:
      tmp[0]=dname[0];
      status=poc_standardaxis(lonname,name, "degree_east","longitude","longitude",grid.xmin,grid.xmax, tmp,cdfgrid->lon);
      status=create_ncvariable(input, cdfgrid->lon);
      tmp[0]=dname[1];
      status=poc_standardaxis(latname,name, "degree_north","latitude","latitude", grid.ymin,grid.ymax, tmp,cdfgrid->lat);
      status=create_ncvariable(input, cdfgrid->lat);
      break;
    default:
      tmp[0]=dname[1];
      tmp[1]=dname[0];
      status=poc_standardaxis_xy(lonname,name, "degree_east","longitude","longitude",grid.xmin,grid.xmax, tmp,cdfgrid->lon);
      status=create_ncvariable(input, cdfgrid->lon);
      status=poc_standardaxis_xy(latname,name, "degree_north","latitude","latitude", grid.ymin,grid.ymax, tmp,cdfgrid->lat);
      status=create_ncvariable(input, cdfgrid->lat);
      if(grid.projection!=0) {
        status=poc_standardaxis_xy(xname,name, "m","x","x", grid.xmin,grid.xmax, tmp,cdfgrid->x);
        status=create_ncvariable(input, cdfgrid->x);
        status=poc_standardaxis_xy(yname,name, "m","y","y", grid.ymin,grid.ymax, tmp,cdfgrid->y);
        status=create_ncvariable(input, cdfgrid->y);
        }
      break;
    }
  
  deletep(&dname[0]);
  deletep(&dname[1]);
  deletep(&lonname);
  deletep(&latname);
  deletep(&levelsname);
  deletep(&xname);
  deletep(&yname);
  
  stat=nc_open(input,NC_WRITE,&ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  if(grid.projection!=0) {
    stat=nc_put_var_double(ncid,cdfgrid->x->id,grid.x);
    nc_check_error(stat,__LINE__,__FILE__);
    stat=nc_put_var_double(ncid,cdfgrid->y->id,grid.y);
    nc_check_error(stat,__LINE__,__FILE__);
    double *lon=new double[grid.Hsize()];
    double *lat=new double[grid.Hsize()];
#pragma omp parallel for
    for (int m=0;m<grid.nx*grid.ny;m++) {
      lat[m]=grid.y[m];
      lon[m]=grid.x[m];
      }
    stat=projection_to_geo(grid.proj4options, lon, lat, grid.Hsize());
    for (int j=0;j<grid.ny;j++) {
      int m=j*grid.nx;
      lon[m]=degree_recale(lon[m],lon[0]);
      for (int i=0;i<grid.nx;i++) {
        int n=j*grid.nx+i;
        lon[n]=degree_recale(lon[n],lon[m]);
        }
      }
    stat=nc_put_var_double(ncid,cdfgrid->lon->id,lon);
    nc_check_error(stat,__LINE__,__FILE__);
    stat=nc_put_var_double(ncid,cdfgrid->lat->id,lat);
    nc_check_error(stat,__LINE__,__FILE__);
    delete[] lon;
    delete[] lat;
    }
  else {
    stat=nc_put_var_double(ncid,cdfgrid->lon->id,grid.x);
    nc_check_error(stat,__LINE__,__FILE__);
    stat=nc_put_var_double(ncid,cdfgrid->lat->id,grid.y);
    nc_check_error(stat,__LINE__,__FILE__);
    }
//   if(grid.z!=NULL) {
//     stat=nc_put_var_double(ncid,cdfgrid->z->id,grid.z);
//     nc_check_error(stat,__LINE__,__FILE__);
//     }

  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  get file information */
  cdfgrid->name     = poc_strdup(name);
  cdfgrid->filename = poc_strdup(input);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_sphericalgrid_xy(const char *input, pocgrd_t *cdfgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  exitIfNull(cdfgrid->x   = new cdfvar_t[1]);
  exitIfNull(cdfgrid->y   = new cdfvar_t[1]);
  exitIfNull(cdfgrid->lon = new cdfvar_t[1]);
  exitIfNull(cdfgrid->lat = new cdfvar_t[1]);
  exitIfNull(cdfgrid->z   = new cdfvar_t[1]);

  status=cdf_varinfo(input, "x",   cdfgrid->x, 0);
  if((status!=0) && (status!=-49)) return(status);
  status=cdf_varinfo(input, "y",   cdfgrid->y, 0);
  if((status!=0) && (status!=-49)) return(status);
  status=cdf_varinfo(input, "z",   cdfgrid->z, 0);
  if((status!=0) && (status!=-49)) return(status);
  status=cdf_varinfo(input, "lon", cdfgrid->lon, 0);
  if((status!=0) && (status!=-49)) return(status);
  status=cdf_varinfo(input, "lat", cdfgrid->lat, 0);
  if((status!=0) && (status!=-49)) return(status);

  cdfgrid->name     = strdup("");
  cdfgrid->filename = strdup(input);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_inq_grid(const char *input, const char *hlabel, const char *vlabel, pocgrd_t *cdfgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  char varname[32];
  
  exitIfNull(cdfgrid->x   = new cdfvar_t[1]);
  exitIfNull(cdfgrid->y   = new cdfvar_t[1]);
  exitIfNull(cdfgrid->lon = new cdfvar_t[1]);
  exitIfNull(cdfgrid->lat = new cdfvar_t[1]);
  exitIfNull(cdfgrid->z   = new cdfvar_t[1]);

  if(strlen(hlabel)!=0) sprintf(varname,"x_%s",hlabel);
  else                  sprintf(varname,"x",   hlabel);
  status=cdf_varinfo(input, varname, cdfgrid->x, 0);
  if((status!=0) && (status!=-49)) return(status);
  
  if(strlen(hlabel)!=0) sprintf(varname,"y_%s",hlabel);
  else                  sprintf(varname,"y",   hlabel);
  status=cdf_varinfo(input, varname, cdfgrid->y, 0);
  if((status!=0) && (status!=-49)) return(status);
  
  if(strlen(vlabel)!=0) sprintf(varname,"depths_%s",vlabel);
  else                   sprintf(varname,"depths",  vlabel);
  status=cdf_varinfo(input, varname, cdfgrid->z, 0);
  if((status!=0) && (status!=-49)) return(status);
  
  if(strlen(hlabel)!=0) sprintf(varname,"lon_%s",hlabel);
  else                  sprintf(varname,"lon",   hlabel);
  status=cdf_varinfo(input, varname, cdfgrid->lon, 0);
  if((status!=0) && (status!=-49)) return(status);
  
  if(strlen(hlabel)!=0) sprintf(varname,"lat_%s",hlabel);
  else                  sprintf(varname,"lat",   hlabel);
  status=cdf_varinfo(input, varname, cdfgrid->lat, 0);
  if((status!=0) && (status!=-49)) return(status);

  cdfgrid->name     = new char[32];
  if( (strlen(hlabel)!=0) &&  (strlen(vlabel)!=0) ) {
    sprintf(cdfgrid->name,"%sx%s",hlabel,vlabel);
    }
  else if(strlen(hlabel)!=0) {
    sprintf(cdfgrid->name,"%s",hlabel);
    }
  else if(strlen(vlabel)!=0) {
    sprintf(cdfgrid->name,"%s",vlabel);
    }
  else {
    sprintf(cdfgrid->name,"%s","");
    }  cdfgrid->filename = strdup(input);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_sphericalgrid_xyt(const char *input, const char *name, grid_t grid, pocgrd_t *cdfgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char     *dname[10];
  char     *tmp[10];
  int dimlgth[10];
  int ncid,stat,dim_id;

  cdfgbl_t global;
  cdfgrid->id=-1;
  char *lonname=NULL, *latname=NULL, *levelsname=NULL;
  bool dim_exist[4];

/*------------------------------------------------------------------------------

  Defines a grid to which variable will be fully or partially connected.

  1) defines dimensions
     beware the unlimited dimension

  2) creates standard horizontal and vertical axis


  3) returns informations to be passed to variable definition and creation routines

*------------------------------------------------------------------------*/

  l=strlen(name);
  if(l==0) {
    dname[0]=strdup("nx");
    dname[1]=strdup("ny");
    dname[2]=strdup("nt");
    lonname=strdup("lon");
    latname=strdup("lat");
    }
  else {
    for(k=0;k<4;k++) {
      exitIfNull(
        dname[k]=(char *) malloc(l+2)
        );
      }
    sprintf(dname[0],"nx_%s",name);
    sprintf(dname[1],"ny_%s",name);
    /*sprintf(dname[3],"t_%s",name);*/
    dname[2]=strdup("nt");
    exitIfNull(
      lonname=(char *) malloc(l+4)
      );
    exitIfNull(
      latname=(char *) malloc(l+4)
      );
    exitIfNull(
      levelsname=(char *) malloc(l+6)
      );
    sprintf(lonname,"lon_%s",name);
    sprintf(latname,"lat_%s",name);
    }

/*------------------------------------------------------------------------------
  enter define mode */
  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  dimlgth[2]=NC_UNLIMITED;
  
  stat = nc_open(input, NC_WRITE, &ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  
  for(k=0;k<2;k++) {
    size_t check;
    stat=nc_inq_dimid(ncid, dname[k], &dim_id);
    if(stat==NC_NOERR) {
      dim_exist[k]=true;
      stat=nc_inq_dim(ncid, dim_id, dname[k], &check);
      if(check!=dimlgth[k]) {
        printf("inconsistent size for dim=%s (%d found, %d wanted)\n",dname[k], check, dimlgth[k]);
        stat = nc_close(ncid);
        return(-1);
        }
      }
    else
      dim_exist[k]=false;
    }
  
  
  stat = nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create dimensions */
  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  dimlgth[2]=NC_UNLIMITED;

  for(k=0;k<2;k++) {
    if(dim_exist[k]) continue;
    stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
    if(stat!=NC_NOERR) {
      printf("nc_def_dim status %d: dimension %s length=%d\n",stat,dname[k], dimlgth[k]);
      nc_check_error(stat,__LINE__,__FILE__);
      }
    }
    
/* *----------------------------------------------------------------------------
  k=3 => dim t :: Check if dim t alreaady created */
  stat=nc_inq_dimid(ncid,dname[k],&dim_id);
  if (stat ==  NC_EBADDIM)  stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
  nc_check_error(stat,__LINE__,__FILE__);

  stat = nc_enddef (ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create axis variable */

  tmp[0]=dname[1];
  tmp[1]=dname[0];

  exitIfNull(cdfgrid->x=new cdfvar_t[1]);
  exitIfNull(cdfgrid->y=new cdfvar_t[1]);

  status=poc_standardaxis_xy(lonname,name, "degree_east","longitude","longitude",-180.0,180.0, tmp,cdfgrid->x);
  status=create_ncvariable(input, cdfgrid->x);
  status=poc_standardaxis_xy(latname,name, "degree_north","latitude","latitude",-90.0,+90.0, tmp,cdfgrid->y);
  status=create_ncvariable(input, cdfgrid->y);

  stat=nc_open(input,NC_WRITE,&ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->x->id,grid.x);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->y->id,grid.y);
  nc_check_error(stat,__LINE__,__FILE__);

  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  get file information */
  status= cdf_globalinfo(input,&global,0);

  cdfgrid->name=strdup(name);
  cdfgrid->filename = strdup(input);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_sphericalgrid_xyzt(const char *input, const char *hlabel, const char *vlabel, grid_t grid, pocgrd_t *cdfgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char     *dname[10];
  char     *tmp[10];
  int dimlgth[10];
  int ncid,stat,dim_id;

  cdfgbl_t global;
  cdfgrid->id=-1;
  char *lonname=NULL, *latname=NULL, *levelsname=NULL;
  bool dim_exist[4];

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  Defines a grid to which variable will be fully or partially connected.

  1) defines dimensions
     beware the unlimited dimension

  2) creates standard horizontal and vertical axis


  3) returns informations to be passed to variable definition and creation routines

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  l=strlen(hlabel);
  if(l==0) {
    dname[0]=strdup("nx");
    dname[1]=strdup("ny");
    lonname=strdup("lon");
    latname=strdup("lat");
    }
  else {
    for(k=0;k<2;k++) {
      exitIfNull(
        dname[k]=(char *) malloc(l+5)
        );
      }
    sprintf(dname[0],"nx_%s",hlabel);
    sprintf(dname[1],"ny_%s",hlabel);
    exitIfNull(
      lonname=(char *) malloc(l+5)
      );
    exitIfNull(
      latname=(char *) malloc(l+5)
      );
    sprintf(lonname,"lon_%s",hlabel);
    sprintf(latname,"lat_%s",hlabel);
    }

  l=strlen(vlabel);
  if(l==0) {
    dname[2]=strdup("nz");
    levelsname=strdup("depths");
    }
  else {
    for(k=2;k<3;k++) {
      exitIfNull(
        dname[k]=(char *) malloc(l+4)
        );
      }
    sprintf(dname[2],"nz_%s",vlabel);
    exitIfNull(
      levelsname=(char *) malloc(l+8)
      );
    sprintf(levelsname,"depths_%s",vlabel);
    }
  dname[3]=strdup("nt");

/*------------------------------------------------------------------------------
  enter define mode */
  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  dimlgth[2]=grid.nz;
  dimlgth[3]=NC_UNLIMITED;
  
  stat = nc_open(input, NC_WRITE, &ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  
  for(k=0;k<3;k++) {
    size_t check;
    stat=nc_inq_dimid(ncid, dname[k], &dim_id);
    if(stat==NC_NOERR) {
      dim_exist[k]=true;
      stat=nc_inq_dim(ncid, dim_id, dname[k], &check);
      if(check!=dimlgth[k]) {
        printf("inconsistent size for dim=%s (%d found, %d wanted)\n",dname[k], check, dimlgth[k]);
        stat = nc_close(ncid);
        return(-1);
        }
      }
    else
      dim_exist[k]=false;
    }
  
  
  stat = nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create dimensions */
  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  dimlgth[2]=grid.nz;
  dimlgth[3]=NC_UNLIMITED;

  for(k=0;k<3;k++) {
    if(dim_exist[k]) continue;
    stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
    if(stat!=NC_NOERR) {
      nc_check_error(stat,__LINE__,__FILE__);
      }
    }
    
/* *----------------------------------------------------------------------------
  k=3 => dim t :: Check if dim t alreaady created */
  stat=nc_inq_dimid(ncid,dname[k],&dim_id);
  if (stat ==  NC_EBADDIM)  stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
  nc_check_error(stat,__LINE__,__FILE__);

  stat = nc_enddef (ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create axis variable */

  tmp[0]=dname[1];
  tmp[1]=dname[0];

  exitIfNull(cdfgrid->x=new cdfvar_t[1]);
  exitIfNull(cdfgrid->y=new cdfvar_t[1]);
  exitIfNull(cdfgrid->z=new cdfvar_t[1]);

  switch (grid.modeH) {
    case 0:
      
    case 1:
      tmp[0]=dname[0];
      status=poc_standardaxis(lonname,hlabel, "degree_east","longitude","longitude",grid.xmin,grid.xmax, tmp,cdfgrid->x);
      status=create_ncvariable(input, cdfgrid->x);
      tmp[0]=dname[1];
      status=poc_standardaxis(latname,hlabel, "degree_north","latitude","latitude", grid.ymin,grid.ymax, tmp,cdfgrid->y);
      status=create_ncvariable(input, cdfgrid->y);
      break;
    default:
      tmp[0]=dname[1];
      tmp[1]=dname[0];
      status=poc_standardaxis_xy(lonname,hlabel, "degree_east","longitude","longitude",grid.xmin,grid.xmax, tmp,cdfgrid->x);
      status=create_ncvariable(input, cdfgrid->x);
      status=poc_standardaxis_xy(latname,hlabel, "degree_north","latitude","latitude", grid.ymin,grid.ymax, tmp,cdfgrid->y);
      status=create_ncvariable(input, cdfgrid->y);
//       if(grid.projection!=0) {
//         status=poc_standardaxis_xy(xname,name, "m","x","x", grid.xmin,grid.xmax, tmp,cdfgrid->x);
//         status=create_ncvariable(input, cdfgrid->x);
//         status=poc_standardaxis_xy(yname,name, "m","y","y", grid.ymin,grid.ymax, tmp,cdfgrid->y);
//         status=create_ncvariable(input, cdfgrid->y);
//         }
      break;
    }
    
//   status=poc_standardaxis_xy(lonname,hlabel, "degree_east","longitude","longitude",-180.0,180.0, tmp,cdfgrid->x);
//   status=create_ncvariable(input, cdfgrid->x);
//   status=poc_standardaxis_xy(latname,hlabel, "degree_north","latitude","latitude",-180.0,180.0, tmp,cdfgrid->y);
//   status=create_ncvariable(input, cdfgrid->y);

  tmp[0]=dname[2];
  tmp[1]=dname[1];
  tmp[2]=dname[0];
  status=poc_standardaxis_z(levelsname, grid.zmask,"m","levels_depths","levels depths",tmp,cdfgrid->z,cdfgrid->x->name,cdfgrid->y->name);
  status=create_ncvariable(input, cdfgrid->z);

  stat=nc_open(input,NC_WRITE,&ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->x->id,grid.x);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->y->id,grid.y);
  nc_check_error(stat,__LINE__,__FILE__);

  if(grid.z!=NULL) stat=nc_put_var_double(ncid,cdfgrid->z->id,grid.z);
  if(stat!=0) {
    nc_check_error(stat,__LINE__,__FILE__);
    }

  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  get file information */
  status= cdf_globalinfo(input,&global,0);

  cdfgrid->name     = new char[32];
  if( (strlen(hlabel)!=0) &&  (strlen(vlabel)!=0) ) {
    if(strcmp(hlabel,vlabel)==0) sprintf(cdfgrid->name,"%s",hlabel);
    else sprintf(cdfgrid->name,"%sx%s",hlabel,vlabel);
    }
  else if(strlen(hlabel)!=0) {
    sprintf(cdfgrid->name,"%s",hlabel);
    }
  else if(strlen(vlabel)!=0) {
    sprintf(cdfgrid->name,"%s",vlabel);
    }
  else {
    sprintf(cdfgrid->name,"%s","");
    }
  cdfgrid->filename = strdup(input);

  cdfgrid->filename = strdup(input);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_sphericalgrid_xyz(const char *input, const char *hlabel, const char *vlabel, grid_t grid, pocgrd_t *cdfgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char     *dname[10];
  char     *tmp[10];
  int dimlgth[10];
  int ncid,stat,dim_id;
  bool dim_exist[3];

  cdfgbl_t global;
  cdfgrid->id=-1;
  char *lonname=NULL, *latname=NULL, *levelsname=NULL;

/*------------------------------------------------------------------------------

  Defines a grid to which variable will be fully or partially connected.

  1) defines dimensions
     beware the unlimited dimension

  2) creates standard horizontal and vertical axis


  3) returns informations to be passed to variable definition and creation routines

*------------------------------------------------------------------------*/

  l=strlen(hlabel);
  if(l==0) {
    dname[0]=strdup("nx");
    dname[1]=strdup("ny");
    lonname=strdup("lon");
    latname=strdup("lat");
    }
  else {
    for(k=0;k<2;k++) {
      exitIfNull(
        dname[k]=(char *) malloc(l+2)
        );
      }
    sprintf(dname[0],"nx_%s",hlabel);
    sprintf(dname[1],"ny_%s",hlabel);
    exitIfNull(
      lonname=(char *) malloc(l+4)
      );
    exitIfNull(
      latname=(char *) malloc(l+4)
      );
    sprintf(lonname,"lon_%s",hlabel);
    sprintf(latname,"lat_%s",hlabel);
    }

  l=strlen(vlabel);
  if(l==0) {
    dname[2]=strdup("nz");
    levelsname=strdup("depths");
    }
  else {
    for(k=2;k<3;k++) {
      exitIfNull(
        dname[k]=(char *) malloc(l+5)
        );
      }
    sprintf(dname[2],"nz_%s",vlabel);
    exitIfNull(
      levelsname=(char *) malloc(l+8)
      );
    sprintf(levelsname,"depths_%s",vlabel);
    }

/*------------------------------------------------------------------------------
  dimension check and creation */
  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  dimlgth[2]=grid.nz;
  
  stat = nc_open(input, NC_WRITE, &ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  
  for(k=0;k<3;k++) {
    size_t check;
    stat=nc_inq_dimid(ncid, dname[k], &dim_id);
    if(stat==NC_NOERR) {
      dim_exist[k]=true;
      stat=nc_inq_dim(ncid, dim_id, dname[k], &check);
      if(check!=dimlgth[k]) {
        printf("inconsistent size for dim=%s (%d found, %d wanted)\n",dname[k], check, dimlgth[k]);
        stat = nc_close(ncid);
        return(-1);
        }
      }
    else
      dim_exist[k]=false;
    }
  
/*------------------------------------------------------------------------------
  enter define mode */
  stat = nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);


/*------------------------------------------------------------------------------
  create dimensions */
  for(k=0;k<3;k++) {
/*if(dimlgth[k]!=0) {*/
    if(dim_exist[k]) continue;
    stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
    if(stat!=NC_NOERR) {
      nc_check_error(stat,__LINE__,__FILE__);
      }
/* }*/
    }

  stat = nc_enddef (ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create axis variable */

  tmp[0]=dname[1];
  tmp[1]=dname[0];

//   exitIfNull(cdfgrid->x=new cdfvar_t[1]);
//   exitIfNull(cdfgrid->y=new cdfvar_t[1]);
  exitIfNull(cdfgrid->lon=new cdfvar_t[1]);
  exitIfNull(cdfgrid->lat=new cdfvar_t[1]);
  exitIfNull(cdfgrid->z=new cdfvar_t[1]);

  switch (grid.modeH) {
    case 0:
      
    case 1:
      tmp[0]=dname[0];
      status=poc_standardaxis(lonname,hlabel, "degree_east","longitude","longitude",grid.xmin,grid.xmax, tmp,cdfgrid->lon);
      status=create_ncvariable(input, cdfgrid->lon);
      tmp[0]=dname[1];
      status=poc_standardaxis(latname,hlabel, "degree_north","latitude","latitude", grid.ymin,grid.ymax, tmp,cdfgrid->lat);
      status=create_ncvariable(input, cdfgrid->lat);
      break;
    default:
      tmp[0]=dname[1];
      tmp[1]=dname[0];
      status=poc_standardaxis_xy(lonname,hlabel, "degree_east","longitude","longitude",grid.xmin,grid.xmax, tmp,cdfgrid->lon);
      status=create_ncvariable(input, cdfgrid->lon);
      status=poc_standardaxis_xy(latname,hlabel, "degree_north","latitude","latitude", grid.ymin,grid.ymax, tmp,cdfgrid->lat);
      status=create_ncvariable(input, cdfgrid->lat);
//       if(grid.projection!=0) {
//         status=poc_standardaxis_xy(xname, hlabel, "m","x","x", grid.xmin,grid.xmax, tmp,cdfgrid->x);
//         status=create_ncvariable(input, cdfgrid->x);
//         status=poc_standardaxis_xy(yname, hlabel, "m","y","y", grid.ymin,grid.ymax, tmp,cdfgrid->y);
//         status=create_ncvariable(input, cdfgrid->y);
//         }
      break;
    }

//   status=poc_standardaxis_xy(lonname,hlabel, "degree_east","longitude","longitude",-180.0,180.0, tmp,cdfgrid->lon);
//   status=create_ncvariable(input, cdfgrid->lon);
//   status=poc_standardaxis_xy(latname,hlabel, "degree_north","latitude","latitude",-180.0,180.0, tmp,cdfgrid->lat);
//   status=create_ncvariable(input, cdfgrid->lat);

  tmp[0]=dname[2];
  tmp[1]=dname[1];
  tmp[2]=dname[0];
  status=poc_standardaxis_z(levelsname, grid.zmask,"m","levels_depths","levels depths",tmp,cdfgrid->z,cdfgrid->lon->name,cdfgrid->lat->name);
  status=create_ncvariable(input, cdfgrid->z);

  stat=nc_open(input,NC_WRITE,&ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->lon->id,grid.x);
  if(stat!=NC_NOERR) {
    nc_check_error(stat,__LINE__,__FILE__);
    }
  
  stat=nc_put_var_double(ncid,cdfgrid->lat->id,grid.y);
  if(stat!=NC_NOERR) {
    nc_check_error(stat,__LINE__,__FILE__);
    }

  if(grid.z!=NULL) {
    stat=nc_put_var_double(ncid,cdfgrid->z->id,grid.z);
    nc_check_error(stat,__LINE__,__FILE__);
    }

  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  get file information */
  status= cdf_globalinfo(input,&global,0);

//  cdfgrid->name=strdup(name);
  cdfgrid->name     = new char[32];
  if( (strlen(hlabel)!=0) &&  (strlen(vlabel)!=0) ) {
    sprintf(cdfgrid->name,"%sx%s",hlabel,vlabel);
    }
  else if(strlen(hlabel)!=0) {
    sprintf(cdfgrid->name,"%s",hlabel);
    }
  else if(strlen(vlabel)!=0) {
    sprintf(cdfgrid->name,"%s",vlabel);
    }
  else {
    sprintf(cdfgrid->name,"%s","");
    }
  cdfgrid->filename = strdup(input);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_sphericalgrid_xyzwt(const char *input, const char *name, int nw,grid_t grid, pocgrd_t *cdfgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char     *dname[10];
  char     *tmp[10];
  int dimlgth[10];
  int ncid,stat,dim_id;

  cdfgbl_t global;
  cdfgrid->id=-1;
  char *lonname=NULL, *latname=NULL, *levelsname=NULL;

/*------------------------------------------------------------------------------
 
  Defines a grid to which variable will be fully or partially connected.

  1) defines dimensions
     beware the unlimited dimension

  2) creates standard horizontal and vertical axis


  3) returns informations to be passed to variable definition and creation routines

*------------------------------------------------------------------------*/

  l=strlen(name);
  if(l==0) {
    dname[0]=strdup("nx");
    dname[1]=strdup("ny");
    dname[2]=strdup("nz");
    dname[3]=strdup("nw");
    dname[4]=strdup("nt");
    lonname=strdup("lon");
    latname=strdup("lat");
    levelsname=strdup("levels");
    }
  else {
    for(k=0;k<4;k++) {
      exitIfNull(
        dname[k]=(char *) malloc(l+2)
        );
      }
    sprintf(dname[0],"nx_%s",name);
    sprintf(dname[1],"ny_%s",name);
    sprintf(dname[2],"nz_%s",name);
    /*sprintf(dname[3],"t_%s",name);*/
    dname[3]=strdup("nw");
    dname[4]=strdup("nt");
    exitIfNull(
      lonname=(char *) malloc(l+4)
      );
    exitIfNull(
      latname=(char *) malloc(l+4)
      );
    exitIfNull(
      levelsname=(char *) malloc(l+6)
      );
    sprintf(lonname,"lon_%s",name);
    sprintf(latname,"lat_%s",name);
    sprintf(levelsname,"depths_%s",name);
    }

/*------------------------------------------------------------------------------
  enter define mode */
  stat = nc_open(input, NC_WRITE, &ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);


/*------------------------------------------------------------------------------
  create dimensions */

  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  dimlgth[2]=grid.nz;
  dimlgth[3]=nw;
  dimlgth[4]=NC_UNLIMITED;

  for(k=0;k<4;k++) {
//    if(dimlgth[k]!=0) {
      stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
      nc_check_error(stat,__LINE__,__FILE__);
//      }
    }
  /*k=3 => dim t :: Check if dim t alreaady created */
  stat=nc_inq_dimid(ncid,dname[k],&dim_id);
  if (stat ==  NC_EBADDIM)  stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
  nc_check_error(stat,__LINE__,__FILE__);

  stat = nc_enddef (ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create axis variable */

  tmp[0]=dname[1];
  tmp[1]=dname[0];

  exitIfNull(cdfgrid->x=new cdfvar_t[1]);
  exitIfNull(cdfgrid->y=new cdfvar_t[1]);
  exitIfNull(cdfgrid->z=new cdfvar_t[1]);
 
  status=poc_standardaxis_xy(lonname,name, "degree_east","longitude","longitude",-180.0,180.0, tmp,cdfgrid->x);
  status=create_ncvariable(input, cdfgrid->x);
  status=poc_standardaxis_xy(latname,name, "degree_north","latitude","latitude",-180.0,180.0, tmp,cdfgrid->y);
  status=create_ncvariable(input, cdfgrid->y);

  tmp[0]=dname[2];
  tmp[1]=dname[1];
  tmp[2]=dname[0];
  status=poc_standardaxis_z(levelsname, grid.zmask,"m","levels_depths","levels depths",  tmp, cdfgrid->z, cdfgrid->x->name, cdfgrid->y->name);
  status=create_ncvariable(input, cdfgrid->z);

  stat=nc_open(input,NC_WRITE,&ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->x->id,grid.x);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->y->id,grid.y);
  nc_check_error(stat,__LINE__,__FILE__);

  if(grid.z!=NULL) stat=nc_put_var_double(ncid,cdfgrid->z->id,grid.z);
  nc_check_error(stat,__LINE__,__FILE__);

  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  get file information */
  status= cdf_globalinfo(input,&global,0);

  cdfgrid->name=strdup(name);
  cdfgrid->filename = strdup(input);

  return(0);
}


/*----------------------------------------------------------------------------*/

int poc_sphericalgrid_xietast(char *input, char *name, grid_t grid, pocgrd_t *cdfgrid)

/*----------------------------------------------------------------------------*/
{

  int k,l,status;
  char     *dname[10];
  char     *tmp[10];
  float mask;
  int dimlgth[10];
  int ncid,stat,dim_id;

  cdfgbl_t global;
  cdfgrid->id=-1;
  char *lonname=NULL, *latname=NULL, *levelsname=NULL;

/*------------------------------------------------------------------------------

  Defines a grid to which variable will be fully or partially connected.

  1) defines dimensions
     beware the unlimited dimension

  2) creates standard horizontal and vertical axis


  3) returns informations to be passed to variable definition and creation routines

*------------------------------------------------------------------------*/

  l=strlen(name);
  if(l==0) {
    dname[0]=strdup("xi");
    dname[1]=strdup("eta");
    dname[2]=strdup("s");
    dname[3]=strdup("x");
    dname[4]=strdup("y");
    dname[5]=strdup("z");
    dname[6]=strdup("time");
    lonname=strdup("lon");
    latname=strdup("lat");
    levelsname=strdup("levels");
    }
  else
    {
    for(k=0;k<7;k++) {
      exitIfNull(
        dname[k]=(char *) malloc(l+4)
        );
    }
    sprintf(dname[0],"xi_%s",name);
    sprintf(dname[1],"eta_%s",name);
    sprintf(dname[2],"s_%s",name);
    sprintf(dname[3],"x_%s",name);
    sprintf(dname[4],"y_%s",name);
    sprintf(dname[5],"z_%s",name);
    /*sprintf(dname[3],"t_%s",name);*/
    dname[6]=strdup("time");
    exitIfNull(
      lonname=(char *) malloc(l+4)
      );
    exitIfNull(
      latname=(char *) malloc(l+4)
      );
    exitIfNull(
      levelsname=(char *) malloc(l+6)
      );
    sprintf(lonname,"lon_%s",name);
    sprintf(latname,"lat_%s",name);
    sprintf(levelsname,"levels_%s",name);
    }

/*------------------------------------------------------------------------------
  enter define mode */
  stat = nc_open(input, NC_WRITE, &ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);


/*------------------------------------------------------------------------------
  create dimensions */

  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  dimlgth[2]=grid.nz;
  dimlgth[3]=grid.nx;
  dimlgth[4]=grid.ny;
  dimlgth[5]=grid.nz;
  dimlgth[6]=NC_UNLIMITED;

  for(k=0;k<6;k++) {
      /*if(dimlgth[k]!=0) {*/
      stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
      nc_check_error(stat,__LINE__,__FILE__);
      /* }*/
    }
  /*k=3 => dim t :: Check if dim t alreaady created */
  stat=nc_inq_dimid(ncid,dname[k],&dim_id);
  if (stat ==  NC_EBADDIM)  stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
  nc_check_error(stat,__LINE__,__FILE__);

  stat = nc_enddef (ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create axis variable */

  tmp[0]=dname[1];
  tmp[1]=dname[0];

  exitIfNull(cdfgrid->x=new cdfvar_t[1]);
  exitIfNull(cdfgrid->y=new cdfvar_t[1]);
  exitIfNull(cdfgrid->z=new cdfvar_t[1]);

  status=poc_standardaxis_xieta(lonname,name, (char *)("degree_east"),(char *)("longitude"),(char *)("longitude"),-180.0,180.0, tmp,cdfgrid->x);
  status=create_ncvariable(input, cdfgrid->x);
  status=poc_standardaxis_xieta(latname,name, (char *)("degree_north"),(char *)("latitude"),(char *)("latitude"),-180.0,180.0, tmp,cdfgrid->y);
  status=create_ncvariable(input, cdfgrid->y);

  tmp[0]=dname[2];
  tmp[1]=dname[1];
  tmp[2]=dname[0];
  mask=1.e+10;
  status=poc_standardaxis_s(levelsname, mask,(char *)("m"),(char *)("levels_depths"),(char *)("levels depths"),
                            tmp,cdfgrid->z,name);
  status=create_ncvariable(input, cdfgrid->z);

  stat=nc_open(input,NC_WRITE,&ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->x->id,grid.x);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->y->id,grid.y);
  nc_check_error(stat,__LINE__,__FILE__);

  if(grid.z!=NULL) stat=nc_put_var_double(ncid,cdfgrid->z->id,grid.z);
  nc_check_error(stat,__LINE__,__FILE__);

  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  get file information */
  status= cdf_globalinfo(input,&global,0);

  cdfgrid->name=strdup(name);
  cdfgrid->filename = strdup(input);

  return(0);
}

/*----------------------------------------------------------------------------*/

int poc_sphericalgrid_xietaswt(char *input, char *name, int nw,grid_t grid, pocgrd_t *cdfgrid)

/*----------------------------------------------------------------------------*/
{

  int k,l,status;
  char     *dname[10];
  char     *tmp[10];
  float mask;
  int dimlgth[10];
  int ncid,stat,dim_id;

  cdfgbl_t global;
  cdfgrid->id=-1;
  char *lonname=NULL, *latname=NULL, *levelsname=NULL;

/*------------------------------------------------------------------------------
 
  Defines a grid to which variable will be fully or partially connected.

  1) defines dimensions
     beware the unlimited dimension

  2) creates standard horizontal and vertical axis


  3) returns informations to be passed to variable definition and creation routines

*------------------------------------------------------------------------*/

  l=strlen(name);
  if(l==0) {
    dname[0]=strdup("xi");
    dname[1]=strdup("eta");
    dname[2]=strdup("s");
    dname[3]=strdup("w");
    dname[4]=strdup("x");
    dname[5]=strdup("y");
    dname[6]=strdup("z");
    dname[7]=strdup("time");
    lonname=strdup("lon");
    latname=strdup("lat");
    levelsname=strdup("levels");
    }
  else
    {
    for(k=0;k<7;k++) {
      exitIfNull(
        dname[k]=(char *) malloc(l+4)
        );
    }
    sprintf(dname[0],"xi_%s",name);
    sprintf(dname[1],"eta_%s",name);
    sprintf(dname[2],"s_%s",name);
    /*sprintf(dname[3],"t_%s",name);*/
    dname[3]=strdup("w");
    sprintf(dname[4],"x_%s",name);
    sprintf(dname[5],"y_%s",name);
    sprintf(dname[6],"z_%s",name);
    dname[7]=strdup("time");
    exitIfNull(
      lonname=(char *) malloc(l+4)
      );
    exitIfNull(
      latname=(char *) malloc(l+4)
      );
    exitIfNull(
      levelsname=(char *) malloc(l+6)
      );
    sprintf(lonname,"lon_%s",name);
    sprintf(latname,"lat_%s",name);
    sprintf(levelsname,"levels_%s",name);
    }

/*------------------------------------------------------------------------------
  enter define mode */
  stat = nc_open(input, NC_WRITE, &ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_redef(ncid);
  nc_check_error(stat,__LINE__,__FILE__);


/*------------------------------------------------------------------------------
  create dimensions */

  dimlgth[0]=grid.nx;
  dimlgth[1]=grid.ny;
  dimlgth[2]=grid.nz;
  dimlgth[3]=nw;
  dimlgth[4]=grid.nx;
  dimlgth[5]=grid.ny;
  dimlgth[6]=grid.nz;
  dimlgth[7]=NC_UNLIMITED;

  for(k=0;k<7;k++) {
      /*if(dimlgth[k]!=0) {*/
      stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
      nc_check_error(stat,__LINE__,__FILE__);
      /* }*/
    }
  /*k=3 => dim t :: Check if dim t alreaady created */
  stat=nc_inq_dimid(ncid,dname[k],&dim_id);
  if (stat ==  NC_EBADDIM)  stat=nc_def_dim(ncid,dname[k], dimlgth[k], &dim_id);
  nc_check_error(stat,__LINE__,__FILE__);

  stat = nc_enddef (ncid);
  nc_check_error(stat,__LINE__,__FILE__);
  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  create axis variable */

  tmp[0]=dname[1];
  tmp[1]=dname[0];

  exitIfNull(cdfgrid->x=new cdfvar_t[1]);
  exitIfNull(cdfgrid->y=new cdfvar_t[1]);
  exitIfNull(cdfgrid->z=new cdfvar_t[1]);
 
  status=poc_standardaxis_xieta(lonname,name, (char *)("degree_east"),(char *)("longitude"),(char *)("longitude"),-180.0,180.0, tmp,cdfgrid->x);
  status=create_ncvariable(input, cdfgrid->x);
  status=poc_standardaxis_xieta(latname,name, (char *)("degree_north"),(char *)("latitude"),(char *)("latitude"),-180.0,180.0, tmp,cdfgrid->y);
  status=create_ncvariable(input, cdfgrid->y);

  tmp[0]=dname[2];
  tmp[1]=dname[1];
  tmp[2]=dname[0];
  mask=1.e+10;
  status=poc_standardaxis_s(levelsname, mask,(char *)("m"),(char *)("levels_depths"),(char *)("levels depths"),
                            tmp,cdfgrid->z,name);
  status=create_ncvariable(input, cdfgrid->z);

  stat=nc_open(input,NC_WRITE,&ncid);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->x->id,grid.x);
  nc_check_error(stat,__LINE__,__FILE__);

  stat=nc_put_var_double(ncid,cdfgrid->y->id,grid.y);
  nc_check_error(stat,__LINE__,__FILE__);

  if(grid.z!=NULL) stat=nc_put_var_double(ncid,cdfgrid->z->id,grid.z);
  nc_check_error(stat,__LINE__,__FILE__);

  stat = nc_close(ncid);
  nc_check_error(stat,__LINE__,__FILE__);

/*------------------------------------------------------------------------------
  get file information */
  status= cdf_globalinfo(input,&global,0);

  cdfgrid->name=strdup(name);
  cdfgrid->filename = strdup(input);

  return(0);
}
