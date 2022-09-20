
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

\brief variable and axis parsing poc-netcdf definitions

Old note : variables definition routines
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"
#include "poc-netcdf.def"
#include "netcdf-proto.h"
#include "geo.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_standardaxis(const char *vname,const char *gname,const char *units,const char *sname,const char *lname,double vmin,double vmax, char **dname, cdfvar_t *variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char *tmpname=NULL;

  variable->id=-1;
  variable->type=NC_DOUBLE;

  variable->initdim(1);


  variable->initatt(5);

  variable->name=strdup(vname);

  k=0;
  variable->dim[k].name=strdup(dname[k]);

  k=0;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("units");
  variable->att[k].data=strdup(units);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("long_name");
  variable->att[k].data=strdup(lname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("standard_name");
  variable->att[k].data=strdup(sname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_DOUBLE;
  variable->att[k].name=strdup("valid_min");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(double))
    );
  for(l=0;l<sizeof(double);l++) variable->att[k].data[l]=(char) *((char *) &vmin +l);
  variable->att[k].length=1;

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_DOUBLE;
  variable->att[k].name=strdup("valid_max");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(double))
    );
  for(l=0;l<sizeof(double);l++) variable->att[k].data[l]=(char) *((char *) &vmax +l);
  variable->att[k].length=1;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_standardaxis_xy(const char *vname,const char *gname,const char *units,const char *sname,const char *lname,double vmin,double vmax, char **dname, cdfvar_t *variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char *tmpname=NULL;

  variable->id=-1;
  variable->type=NC_DOUBLE;

  variable->initdim(2);


  variable->initatt(7);

  variable->name=poc_strdup(vname);

  k=0;
  variable->dim[k].name=poc_strdup(dname[k]);
  k++;
  variable->dim[k].name=poc_strdup(dname[k]);

  k=0;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=poc_strdup("units");
  variable->att[k].data=poc_strdup(units);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=poc_strdup("long_name");
  variable->att[k].data=poc_strdup(lname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=poc_strdup("standard_name");
  variable->att[k].data=poc_strdup(sname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable->att[k].name=poc_strdup("axis");
  variable->att[k].name=poc_strdup("content");
  variable->att[k].data=poc_strdup("YX");
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=poc_strdup("associate");
  l = strlen(gname);
  if (l >0) {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 10)
      );
    sprintf(tmpname,"lat_%s lon_%s", gname,gname);
    }
  else {
    tmpname=strdup("lat lon");
    }
  variable->att[k].data=poc_strdup(tmpname);
  variable->att[k].length=strlen( variable->att[k].data);
  free(tmpname);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_DOUBLE;
  variable->att[k].name=poc_strdup("valid_min");
  exitIfNull(
    variable->att[k].data=new char[sizeof(double)]
    );
  for(l=0;l<sizeof(double);l++) variable->att[k].data[l]=(char) *((char *) &vmin +l);
  variable->att[k].length=1;

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_DOUBLE;
  variable->att[k].name=poc_strdup("valid_max");
  exitIfNull(
    variable->att[k].data=new char[sizeof(double)]
    );
  for(l=0;l<sizeof(double);l++) variable->att[k].data[l]=(char) *((char *) &vmax +l);
  variable->att[k].length=1;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_standardaxis_z(const char *name,float mask,const char *units,const char *standardname,const char *longname, char **dname, cdfvar_t *variable, const char *xname, const char *yname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  float scale=1.0, offset=0.0;
  int k,l,status;
  char *tmp=NULL;

  variable->id=-1;
  variable->type=NC_FLOAT;

  variable->initdim(3);

  variable->initatt(10);

  variable->name=strdup(name);

  k=0;
  variable->dim[k].name=strdup(dname[k]);
  k++;
  variable->dim[k].name=strdup(dname[k]);
  k++;
  variable->dim[k].name=strdup(dname[k]);

  k=0;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("units");
  variable->att[k].data=strdup(units);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("long_name");
  variable->att[k].data=strdup(longname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("standard_name");
  variable->att[k].data=strdup(standardname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("short_name");
  variable->att[k].data=strdup(name);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable->att[k].name=strdup("axis");
  variable->att[k].name=strdup("content");
  variable->att[k].data=strdup("ZYX");
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("associate");
//   l=strlen(gridname);
//   if(l==0) {
//     variable->att[k].data=strdup("level lat lon");
//     }
//   else {
//     if ((tmp = (char *) malloc((l*3+100)*sizeof(char)))==NULL) {
//       __ERR_BASE_LINE__("");perror("tmp");
//       exit(-1);
//       }
//     sprintf(tmp,"levels_%s lat_%s lon_%s ", gridname,gridname,gridname);
//     variable->att[k].data=strdup(tmp);
//     free(tmp);
//     }
  exitIfNull(
    tmp = (char *) malloc(strlen(name)+strlen(yname)+strlen(xname)+3)
    );
  sprintf(tmp,"%s %s %s", name,yname,xname);
  variable->att[k].data=strdup(tmp);
  variable->att[k].length=strlen( variable->att[k].data);
  free(tmp);

  k++;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=strdup("missing_value");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &mask +l);
  variable->att[k].length=1;

  k++;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=strdup("_FillValue");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &mask +l);
  variable->att[k].length=1;

  k++;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=strdup("scale_factor");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable->att[k].data,&scale,sizeof(float));
  variable->att[k].length=1;

  k++;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=strdup("add_offset"); 
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &offset +l);
  variable->att[k].length=1;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poc_standardaxis_incidence(const char *name, const char *standardname,const char *longname, const char **dname, cdfvar_t *variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  float scale=1.0, offset=0.0;
  int k,l,status;

  variable->id=-1;
  variable->type=NC_INT;

  variable->initdim(1);
  variable->initatt(5);

  variable->name=strdup(name);

  k=0;
  variable->dim[k].name=strdup(dname[k]);

  k=0;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("long_name");
  variable->att[k].data=strdup(longname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("standard_name");
  variable->att[k].data=strdup(standardname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("short_name");
  variable->att[k].data=strdup(name);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable->att[k].name=strdup("axis");
  variable->att[k].name=strdup("content");
  variable->att[k].data=strdup("ZYX");
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("associate");
  variable->att[k].data=strdup("level lat lon");
  variable->att[k].length=strlen( variable->att[k].data);

   return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

cdfvar_t poc_standardtime(const char *vname,const char *dname,const char *units,date_t origin)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  cdfvar_t variable;
  int k,l,status;
  char *s=NULL;
  char tmp[1024];

  variable.id=-1;
  variable.type=NC_DOUBLE;

  variable.initdim(1);

  variable.initatt(6);

  variable.name=strdup(vname);

  variable.dim[0].name=strdup(dname);

  k=0;
  s=poc_getdate(origin);
  sprintf(tmp,"%s from %s",units,s);
  
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(tmp);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("time_origin");
  variable.att[k].data=strdup(s);
  variable.att[k].length=strlen( variable.att[k].data);
  free(s);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("calendar");
  variable.att[k].data=strdup("gregorian");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup("time");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("T");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  variable.att[k].data=strdup("time");
  variable.att[k].length=strlen( variable.att[k].data);

  return(variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void poc_standardvariable_xy(cdfvar_t *variable,const char *name, float mask,const char *units,float scale,float offset, const char *standardname,const char *longname,const char *shortname,pocgrd_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char *tmpname=NULL;

  variable->id=-1;
  variable->type=NC_FLOAT;

  variable->initdim(2);

  variable->initatt(10);

  variable->name=poc_strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable->dim[k].name=poc_strdup("ny");
    k++;
    variable->dim[k].name=poc_strdup("nx");
    }
  else {
    exitIfNull( tmpname=new char[l+3] );
    k=0;
    sprintf(tmpname,"ny_%s",grid.name);
    variable->dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"nx_%s",grid.name);
    variable->dim[k].name=strdup(tmpname);
    free(tmpname);
    }

  k=0;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=poc_strdup("units");
  variable->att[k].data=poc_strdup(units);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=poc_strdup("long_name");
  variable->att[k].data=poc_strdup(longname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=poc_strdup("standard_name");
  variable->att[k].data=poc_strdup(standardname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=poc_strdup("short_name");
  variable->att[k].data=poc_strdup(shortname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable->att[k].name=strdup("axis");
  variable->att[k].name=poc_strdup("content");
  variable->att[k].data=poc_strdup("YX");
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=poc_strdup("associate");
  if(l==0) {
    tmpname =  poc_strdup("lat lon");
    }
  else {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 10)
      );
    sprintf(tmpname,"lat_%s lon_%s", grid.name,grid.name);
    }
  variable->att[k].data=poc_strdup(tmpname);
  variable->att[k].length=strlen( variable->att[k].data);
  delete[] tmpname;

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=poc_strdup("missing_value");
  exitIfNull(
    variable->att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &mask +l);
  variable->att[k].length=1;

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=poc_strdup("_FillValue");
  exitIfNull(
    variable->att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &mask +l);
  variable->att[k].length=1;

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=poc_strdup("scale_factor");
  exitIfNull(
    variable->att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable->att[k].data,&scale,sizeof(float));
  variable->att[k].length=1;

  k++;
  variable->att[k].id=k;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=poc_strdup("add_offset");
  exitIfNull(
    variable->att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &offset +l);
  variable->att[k].length=1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_standardvariable_xyt(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname,pocgrd_t grid, cdfvar_t &variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,l,status;
  char *tmpname;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(3);

  variable.initatt(10);

  variable.name=new char[strlen(name)+1];
  strcpy(variable.name,name);

//  variable.name=poc_strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=poc_strdup("nt");
    k++;
    variable.dim[k].name=poc_strdup("ny");
    k++;
    variable.dim[k].name=poc_strdup("nx");
    }
  else {
    tmpname = (char *) malloc(l+3);
    k=0;

    sprintf(tmpname,"nt");
    variable.dim[k].name=poc_strdup(tmpname);
    k++;
    sprintf(tmpname,"ny_%s",grid.name);
    variable.dim[k].name=poc_strdup(tmpname);

    k++;
    sprintf(tmpname,"nx_%s",grid.name);
    variable.dim[k].name=poc_strdup(tmpname);
    free(tmpname);
    }

  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("units");
  variable.att[k].data=poc_strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("long_name");
  variable.att[k].data=poc_strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("standard_name");
  variable.att[k].data=poc_strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("short_name");
  variable.att[k].data=poc_strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("content");
  variable.att[k].data=poc_strdup("TYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("associate");
  if(l==0) {
     tmpname =  strdup("time lat lon");
    }
  else {
    tmpname = (char *) malloc(2*l + 15);
    sprintf(tmpname,"time lat_%s lon_%s",grid.name,grid.name);
    }
  variable.att[k].data=poc_strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("missing_value");
  variable.att[k].data=new char[sizeof(float)];
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("_FillValue");
  variable.att[k].data=new char[sizeof(float)];
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("scale_factor");
  variable.att[k].data=new char[sizeof(float)];
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data,&scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("add_offset");
  variable.att[k].data=new char[sizeof(float)];
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int poc_standardvariable_xyt(const char *name, double mask,const char *units,double scale,double offset,
			      const char *standardname,const char *longname,const char *shortname,pocgrd_t grid, cdfvar_t &variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   cdfvar_t variable;
  int k,l,status;
  char *tmpname;

  variable.id=-1;
  variable.type=NC_DOUBLE;

  variable.initdim(3);

  variable.initatt(10);

  variable.name=new char[strlen(name)+1];
  strcpy(variable.name,name);

//  variable.name=poc_strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("nt");
    k++;
    variable.dim[k].name=strdup("ny");
    k++;
    variable.dim[k].name=strdup("nx");
    }
  else {
    tmpname = (char *) malloc(l+3);
    k=0;

    sprintf(tmpname,"nt");
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"ny_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);

    k++;
    sprintf(tmpname,"nx_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
    }

  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("time lat lon");
    }
  else {
    tmpname = (char *) malloc(2*l + 15);
    sprintf(tmpname,"time lat_%s lon_%s",grid.name,grid.name);
    }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_DOUBLE;
  variable.att[k].name=strdup("missing_value");
  variable.att[k].data=new char[sizeof(double)];
  for(l=0;l<sizeof(double);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_DOUBLE;
  variable.att[k].name=strdup("_FillValue");
  variable.att[k].data=new char[sizeof(double)];
  for(l=0;l<sizeof(double);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_DOUBLE;
  variable.att[k].name=strdup("scale_factor");
  variable.att[k].data=new char[sizeof(double)];
  for(l=0;l<sizeof(double);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data,&scale,sizeof(double));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_DOUBLE;
  variable.att[k].name=strdup("add_offset");
  variable.att[k].data=new char[sizeof(double)];
  for(l=0;l<sizeof(double);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_standardvariable_xyt(const char *name, short mask,const char *units,float scale,float offset,
                               const char *standardname, const char *longname, const char *shortname,pocgrd_t grid, cdfvar_t &variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   cdfvar_t variable;
  int k,l,status;
  char *tmpname;

  variable.id=-1;
  variable.type=NC_SHORT;

  variable.initdim(3);

  variable.initatt(10);

  variable.name=new char[strlen(name)+1];
  strcpy(variable.name,name);

//  variable.name=poc_strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("nt");
    k++;
    variable.dim[k].name=strdup("ny");
    k++;
    variable.dim[k].name=strdup("nx");
    }
  else {
    tmpname = (char *) malloc(l+3);
    k=0;

    sprintf(tmpname,"nt");
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"ny_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);

    k++;
    sprintf(tmpname,"nx_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
    }

  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("time lat lon");
    }
  else {
    tmpname = (char *) malloc(2*l + 15);
    sprintf(tmpname,"time lat_%s lon_%s",grid.name,grid.name);
    }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_SHORT;
  variable.att[k].name=strdup("missing_value");
  variable.att[k].data=new char[sizeof(short)];
  for(l=0;l<sizeof(short);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_SHORT;
  variable.att[k].name=strdup("_FillValue");
  variable.att[k].data=new char[sizeof(short)];
  for(l=0;l<sizeof(short);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  variable.att[k].data=new char[sizeof(float)];
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data,&scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  variable.att[k].data=new char[sizeof(float)];
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

// OBSOLETE (replaced by newer version from tugo )
// cdfvar_t poc_standardvariable_xyt(const char *name, float mask,const char *units,float scale,float offset,
// 				  const char *standardname,const char *longname,const char *shortname,pocgrd_t grid)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   cdfvar_t variable;
//   int k,l,status;
//   char *tmpname=NULL;
//
//   status=init_ncvariable(&variable);
//
//   variable.id=-1;
//   variable.type=NC_FLOAT;
//
//   variable.initdim(3);
//   variable.dimids =new int[variable.ndim];
//   variable.dimname=new char *[variable.ndim];;
//
//   variable.initatt(10);
//
//   variable.name=strdup(name);
//
//   l = strlen(grid.name);
//   if(l==0) {
//     k=0;
//     variable.dim[k].name=strdup("t");
//     k++;
//     variable.dim[k].name=strdup("y");
//     k++;
//     variable.dim[k].name=strdup("x");
//   }
//   else {
//     if ((tmpname = (char *) malloc(l+3))==NULL) {
//       perror("tmpname");
//       __ERR_BASE_LINE__("exiting\n");exit(-1);
//     }
//     k=0;
//
//     sprintf(tmpname,"t");
//     variable.dim[k].name=strdup(tmpname);
//     k++;
//     sprintf(tmpname,"y_%s",grid.name);
//     variable.dim[k].name=strdup(tmpname);
//
//     k++;
//     sprintf(tmpname,"x_%s",grid.name);
//     variable.dim[k].name=strdup(tmpname);
//     free(tmpname);
//   }
//
//   k=0;
//   variable.att[k].type=NC_CHAR;
//   variable.att[k].name=strdup("units");
//   variable.att[k].data=strdup(units);
//   variable.att[k].length=strlen( variable.att[k].data);
//
//   k++;
//   variable.att[k].type=NC_CHAR;
//   variable.att[k].name=strdup("long_name");
//   variable.att[k].data=strdup(longname);
//   variable.att[k].length=strlen( variable.att[k].data);
//
//   k++;
//   variable.att[k].type=NC_CHAR;
//   variable.att[k].name=strdup("standard_name");
//   variable.att[k].data=strdup(standardname);
//   variable.att[k].length=strlen( variable.att[k].data);
//
//   k++;
//   variable.att[k].type=NC_CHAR;
//   variable.att[k].name=strdup("short_name");
//   variable.att[k].data=strdup(shortname);
//   variable.att[k].length=strlen( variable.att[k].data);
//
//   k++;
//   variable.att[k].type=NC_CHAR;
// /* *------------------------------------------------------------------------
//   attribute name changed to comply with CF standard*/
// //  variable.att[k].name=strdup("axis");
//   variable.att[k].name=strdup("content");
//   variable.att[k].data=strdup("TYX");
//   variable.att[k].length=strlen( variable.att[k].data);
//
//   k++;
//   variable.att[k].type=NC_CHAR;
//   variable.att[k].name=strdup("associate");
//   if(l==0) {
//      tmpname =  strdup("time levels lat lon");
//   }
//   else {
//     if ((tmpname = (char *) malloc(2*l + 15))==NULL) {
//       perror("tmpname");
//       __ERR_BASE_LINE__("exiting\n");exit(-1);
//     }
//     sprintf(tmpname,"time lat_%s lon_%s",
// 	    grid.name,grid.name);
//   }
//   variable.att[k].data=strdup(tmpname);
//   variable.att[k].length=strlen( variable.att[k].data);
//   free(tmpname);
//
//   k++;
//   variable.att[k].type=NC_FLOAT;
//   variable.att[k].name=strdup("missing_value");
//   if ((variable.att[k].data=(char *) malloc(sizeof(float)))==NULL) {
//     perror("variable.att[k].data");
//     __ERR_BASE_LINE__("exiting\n");exit(-1);
//   }
//   for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
//   variable.att[k].length=1;
//
//   k++;
//   variable.att[k].type=NC_FLOAT;
//   variable.att[k].name=strdup("_FillValue");
//   if ((variable.att[k].data=(char *) malloc(sizeof(float)))==NULL) {
//     perror("variable.att[k].data");
//     __ERR_BASE_LINE__("exiting\n");exit(-1);
//   }
//   for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
//   variable.att[k].length=1;
//
//   k++;
//   variable.att[k].type=NC_FLOAT;
//   variable.att[k].name=strdup("scale_factor");
//   if ((variable.att[k].data=(char *) malloc(sizeof(float)))==NULL) {
//     perror("variable.att[k].data");
//     __ERR_BASE_LINE__("exiting\n");exit(-1);
//   }
//   for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
//   memcpy(variable.att[k].data,&scale,sizeof(float));
//   variable.att[k].length=1;
//
//   k++;
//   variable.att[k].type=NC_FLOAT;
//   variable.att[k].name=strdup("add_offset");
//   if ((variable.att[k].data=(char *) malloc(sizeof(float)))==NULL) {
//     perror("variable.att[k].data");
//     __ERR_BASE_LINE__("exiting\n");exit(-1);
//   }
//   for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
//   variable.att[k].length=1;
//
//    return(variable);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

cdfvar_t poc_standardvariable_xyw(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname,pocgrd_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  cdfvar_t variable;
  int k,l,status;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(3);

  variable.initatt(10);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("w");
    k++;
    variable.dim[k].name=strdup("y");
    k++;
    variable.dim[k].name=strdup("x");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+2)
      );
    k=0;
    
    sprintf(tmpname,"w");
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"y_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"x_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  


  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("WYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("freq lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 15)
      );
    sprintf(tmpname,"freq lat_%s lon_%s",
            grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data,&scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_standardvariable_xyz(const char *name, float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname,pocgrd_t grid, cdfvar_t & variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  /*cdfvar_t variable*/;
  int k,l,status;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(3);

  variable.initatt(11);

  variable.name=poc_strdup(name);

//   l = strlen(grid.name);
//   if(l==0) {
//     k=0;
//     variable.dim[k].name=poc_strdup("z");
//     k++;
//     variable.dim[k].name=poc_strdup("y");
//     k++;
//     variable.dim[k].name=poc_strdup("x");
//     }
//   else {
//     if ((tmpname = (char *) malloc(l+3))==NULL) {
//       __ERR_BASE_LINE__("");perror("tmpname");
//       exit(-1);
//       }
//     k=0;
//     sprintf(tmpname,"z_%s",grid.name);
//     variable.dim[k].name=poc_strdup(tmpname);
//     k++;
//     sprintf(tmpname,"y_%s",grid.name);
//     variable.dim[k].name=poc_strdup(tmpname);
//
//     k++;
//     sprintf(tmpname,"x_%s",grid.name);
//     variable.dim[k].name=poc_strdup(tmpname);
//     free(tmpname);
//     }
    
  k=0;
  variable.dim[k].name=poc_strdup(grid.z->dim[0].name);
  k++;
  variable.dim[k].name=poc_strdup(grid.z->dim[1].name);
  k++;
  variable.dim[k].name=poc_strdup(grid.z->dim[2].name);
  
  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("units");
  variable.att[k].data=poc_strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("long_name");
  variable.att[k].data=poc_strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("standard_name");
  variable.att[k].data=poc_strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("short_name");
  variable.att[k].data=poc_strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=poc_strdup("axis");
  variable.att[k].name=poc_strdup("content");
  variable.att[k].data=poc_strdup("ZYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("associate");
  exitIfNull(
    tmpname = (char *) malloc(strlen(grid.z->name)+strlen(grid.lat->name)+strlen(grid.lon->name)+3)
    );
  sprintf(tmpname,"%s %s %s", grid.z->name,grid.lat->name,grid.lon->name);
  variable.att[k].data=poc_strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  
  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("coordinates");
  variable.att[k].data=poc_strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("missing_value");
  exitIfNull(
    variable.att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data,&scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("add_offset");
  exitIfNull(
    variable.att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poc_standardvariable_xyzt(const char *name, float mask,const char *units,float scale,float offset,
                                   const char *standardname,const char *longname,const char *shortname,pocgrd_t grid,cdfvar_t & variable)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  cdfvar_t v1;
  int k,l;
  char *tmpname=NULL;
  
//   variable=new cdfvar_t;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(4);

  variable.initatt(10);
  variable.name=poc_strdup(name);

  l = strlen(grid.name);
//   if(l==0) {
//     k=0;
//     variable.dim[k].name=poc_strdup("nt");
//     k++;
//     variable.dim[k].name=poc_strdup("nz");
//     k++;
//     variable.dim[k].name=poc_strdup("ny");
//     k++;
//     variable.dim[k].name=poc_strdup("nx");
//     }
//   else {
//     exitIfNull(
//       tmpname=(char *) malloc(l+3)
//       );
//     k=0;
//     
//     sprintf(tmpname,"nt");
//     variable.dim[k].name=poc_strdup(tmpname);
//     k++;
//                         
//     sprintf(tmpname,"nz_%s",grid.name);
//     variable.dim[k].name=poc_strdup(tmpname);
//     k++;
//     sprintf(tmpname,"ny_%s",grid.name);
//     variable.dim[k].name=poc_strdup(tmpname);
//     
//     k++;
//     sprintf(tmpname,"nx_%s",grid.name);
//     variable.dim[k].name=poc_strdup(tmpname);
//     free(tmpname);
//     }
  
  k=0;
  variable.dim[k].name=poc_strdup("nt");
  k++;
  variable.dim[k].name=poc_strdup(grid.z->dim[0].name);
  k++;
  variable.dim[k].name=poc_strdup(grid.z->dim[1].name);
  k++;
  variable.dim[k].name=poc_strdup(grid.z->dim[2].name);
  
  
  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("units");
  variable.att[k].data=poc_strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("long_name");
  variable.att[k].data=poc_strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("standard_name");
  variable.att[k].data=poc_strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("short_name");
  variable.att[k].data=poc_strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=poc_strdup("axis");
  variable.att[k].name=poc_strdup("content");
  variable.att[k].data=poc_strdup("TZYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=poc_strdup("associate");
  if(l==0) {
    tmpname =  strdup("time depths lat lon");
    }
  else {
    exitIfNull(
      tmpname=(char *) malloc(3*l + 23)
      );
    sprintf(tmpname,"time depths_%s lat_%s lon_%s", grid.name,grid.name,grid.name);
    }
  variable.att[k].data=poc_strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("missing_value");
  exitIfNull(
    variable.att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=poc_strdup("add_offset");
  exitIfNull(
    variable.att[k].data=new char[sizeof(float)]
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

cdfvar_t poc_standardvariable_xywt(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname,pocgrd_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  cdfvar_t variable,v1;
  int k,l;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(4);

  variable.initatt(10);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("t");
    k++;
    variable.dim[k].name=strdup("w");
    k++;
    variable.dim[k].name=strdup("y");
    k++;
    variable.dim[k].name=strdup("x");
    }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+3)
      );
    k=0;
    
    sprintf(tmpname,"t");
    variable.dim[k].name=strdup(tmpname);
    k++;
                        
    sprintf(tmpname,"w");
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"y_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"x_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
    }
  
  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TWYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
    tmpname =  strdup("time freq lat lon");
    }
  else {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 19)
      );
    sprintf(tmpname,"time freq lat_%s lon_%s",grid.name,grid.name,grid.name);
    }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

  return(variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

cdfvar_t poc_standardvariable_xyzw(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname,pocgrd_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  cdfvar_t variable,v1;
  int k,l;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(4);

  variable.initatt(10);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("w");
    k++;
    variable.dim[k].name=strdup("z");
    k++;
    variable.dim[k].name=strdup("y");
    k++;
    variable.dim[k].name=strdup("x");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+3)
      );
    k=0;
    sprintf(tmpname,"w");
    variable.dim[k].name=strdup(tmpname);

    k++;
    sprintf(tmpname,"z_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
 
    k++;
    sprintf(tmpname,"y_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"x_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  
  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TWYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("time level lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 20)
      );
    sprintf(tmpname,"time level lat_%s lon_%s",
            grid.name,grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

cdfvar_t poc_standardvariable_xyzwt(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname,pocgrd_t grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  cdfvar_t variable,v1;
  int k,l;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(5);

  variable.initatt(10);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("t");
    k++;
    variable.dim[k].name=strdup("w");
    k++;
    variable.dim[k].name=strdup("z");
    k++;
    variable.dim[k].name=strdup("y");
    k++;
    variable.dim[k].name=strdup("x");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+3)
      );
    k=0;
    
    sprintf(tmpname,"t");
    variable.dim[k].name=strdup(tmpname);
    k++;
                        
    sprintf(tmpname,"w");
    variable.dim[k].name=strdup(tmpname);
    k++;

    sprintf(tmpname,"z_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"y_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"x_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  
  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TWZYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("time freq levels lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(3*l + 27)
      );
    sprintf(tmpname,"time freq levels_%s lat_%s lon_%s",
            grid.name,grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

cdfvar_t poc_standardvariable_ift(const char *name,float mask,const char *units,float scale,float offset,const char *standardname,const char *longname,const char *shortname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  cdfvar_t variable;
  int k,l;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(4);

  variable.initatt(10);

  variable.name=strdup(name);

  k=0;
  variable.dim[k].name=strdup("t");
  k++;
  variable.dim[k].name=strdup("z");
  k++;
  variable.dim[k].name=strdup("y");
  k++;
  variable.dim[k].name=strdup("x");

  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TFYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  variable.att[k].data=strdup("time frequency incidence");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}

/*----------------------------------------------------------------------------*/

int poc_standardaxis_xieta(char *vname,char *gname, char *units,char *sname,char *lname,double vmin,double vmax, char **dname, cdfvar_t *variable)

/*----------------------------------------------------------------------------*/
{

  int k,l,status;
  char *tmpname=NULL;

  variable->id=-1;
  variable->type=NC_DOUBLE;

  variable->initdim(2);

  variable->initatt(8);

  variable->name=strdup(vname);

  k=0;
  variable->dim[k].name=strdup(dname[k]);
  k++;
  variable->dim[k].name=strdup(dname[k]);

  k=0;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("units");
  variable->att[k].data=strdup(units);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("long_name");
  variable->att[k].data=strdup(lname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("standard_name");
  variable->att[k].data=strdup(sname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("axis-roms");
  variable->att[k].data=strdup("ETAXI");
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable->att[k].name=strdup("axis");
  variable->att[k].name=strdup("content");
  variable->att[k].data=strdup("YX");
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("associate");
  l = strlen(gname);
  if (l >0) {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 10)
      );
    sprintf(tmpname,"lat_%s lon_%s",
            gname,gname);
  }
  else {
    tmpname=strdup("lat lon");
  }
  variable->att[k].data=strdup(tmpname);
  variable->att[k].length=strlen( variable->att[k].data);
  free(tmpname);

  k++;
  variable->att[k].type=NC_DOUBLE;
  variable->att[k].name=strdup("valid_min");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(double))
    );
  for(l=0;l<sizeof(double);l++) variable->att[k].data[l]=(char) *((char *) &vmin +l);
  variable->att[k].length=1;

  k++;
  variable->att[k].type=NC_DOUBLE;
  variable->att[k].name=strdup("valid_max");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(double))
    );
  for(l=0;l<sizeof(double);l++) variable->att[k].data[l]=(char) *((char *) &vmax +l);
  variable->att[k].length=1;

  return(0);
}

/*----------------------------------------------------------------------------*/

int poc_standardaxis_s(char *name,float mask,char *units,char *standardname,char *longname, char **dname, cdfvar_t *variable, char *gridname)

/*----------------------------------------------------------------------------*/
{

  float scale=1.0, offset=0.0;
  int k,l,status;
  char *tmp=NULL;

  variable->id=-1;
  variable->type=NC_FLOAT;

  variable->initdim(3);
 
  variable->initatt(11);

  variable->name=strdup(name);

  k=0;
  variable->dim[k].name=strdup(dname[k]);
  k++;
  variable->dim[k].name=strdup(dname[k]);
  k++;
  variable->dim[k].name=strdup(dname[k]);

  k=0;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("units");
  variable->att[k].data=strdup(units);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("long_name");
  variable->att[k].data=strdup(longname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("standard_name");
  variable->att[k].data=strdup(standardname);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("short_name");
  variable->att[k].data=strdup(name);
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("axis-roms");
  variable->att[k].data=strdup("SETAXI");
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable->att[k].name=strdup("axis");
  variable->att[k].name=strdup("content");
  variable->att[k].data=strdup("ZYX");
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_CHAR;
  variable->att[k].name=strdup("associate");
  l=strlen(gridname);
  if(l==0) {
      variable->att[k].data=strdup("level lat lon");
    }
  else
    {
      exitIfNull(
        tmp=(char *) malloc((l*3+100)*sizeof(char))
        );
      sprintf(tmp,"levels_%s lat_%s lon_%s ",
        gridname,gridname,gridname);
      variable->att[k].data=strdup(tmp);
      free(tmp);
    }
  variable->att[k].length=strlen( variable->att[k].data);

  k++;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=strdup("missing_value");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &mask +l);
  variable->att[k].length=1;

  k++;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=strdup("_FillValue");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &mask +l);
  variable->att[k].length=1;

  k++;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=strdup("scale_factor");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable->att[k].data,&scale,sizeof(float));
  variable->att[k].length=1;

  k++;
  variable->att[k].type=NC_FLOAT;
  variable->att[k].name=strdup("add_offset");
  exitIfNull(
    variable->att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable->att[k].data[l]=(char) *((char *) &offset +l);
  variable->att[k].length=1;

  return(0);
}



/*----------------------------------------------------------------------------*/

cdfvar_t poc_standardvariable_xieta(char *name, float mask,char *units,float scale,float offset,
                                  char *standardname,char *longname,char *shortname,pocgrd_t grid)

/*----------------------------------------------------------------------------*/
{

  cdfvar_t variable;
  int k,l,status;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(2);

  variable.initatt(11);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("eta");
    k++;
    variable.dim[k].name=strdup("xi");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+3)
      );
    k=0;
    sprintf(tmpname,"eta_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"xi_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  


  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("axis-roms");
  variable.att[k].data=strdup("ETAXI");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("YX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 10)
      );
    sprintf(tmpname,"lat_%s lon_%s",
            grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data,&scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}

/*----------------------------------------------------------------------------*/

cdfvar_t poc_standardvariable_xietat(char *name, float mask,char *units,float scale,float offset,
                                  char *standardname,char *longname,char *shortname,pocgrd_t grid)

/*----------------------------------------------------------------------------*/
{

  cdfvar_t variable;
  int k,l,status;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(3);

  variable.initatt(11);
  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("time");
    k++;
    variable.dim[k].name=strdup("eta");
    k++;
    variable.dim[k].name=strdup("xi");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+5)
      );
    k=0;
    
    sprintf(tmpname,"time");
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"eta_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"xi_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  
  /*printf("   -1-*** std var variable->name=%s \n",variable.name);*/
  /*for (k=0;k<variable.ndim;k++)*/
    /*printf("   std var variable->dimname[%d]=%s \n",k,variable.dim[k].name);*/



  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("axis-roms");
  variable.att[k].data=strdup("TETAXI");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("time levels lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 15)
      );
    sprintf(tmpname,"time lat_%s lon_%s",
            grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data,&scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;
  /*printf("   *** std var variable->name=%s \n",variable.name);*/
  /*for (k=0;k<variable.ndim;k++)
    printf("   std var variable->dimname[%d]=%s \n",k,variable.dim[k].name);*/

   return(variable);
}


/*----------------------------------------------------------------------------*/

cdfvar_t poc_standardvariable_xietaw(char *name, float mask,char *units,float scale,float offset,
                                  char *standardname,char *longname,char *shortname,pocgrd_t grid)

/*----------------------------------------------------------------------------*/
{

  cdfvar_t variable;
  int k,l,status;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(3);

  variable.initatt(11);
  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("w");
    k++;
    variable.dim[k].name=strdup("eta");
    k++;
    variable.dim[k].name=strdup("xi");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+2)
      );
    k=0;
    
    sprintf(tmpname,"w");
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"eta_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"xi_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  


  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("axis-roms");
  variable.att[k].data=strdup("WETAX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("WYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("freq lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 15)
      );
    sprintf(tmpname,"freq lat_%s lon_%s",
            grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data,&scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}


/*----------------------------------------------------------------------------*/

cdfvar_t poc_standardvariable_xietas(char *name, float mask,char *units,float scale,float
                                  offset,char *standardname,char *longname,
                                  char *shortname,pocgrd_t grid)

/*----------------------------------------------------------------------------*/
{

  cdfvar_t variable;
  int k,l,status;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(3);

  variable.initatt(11);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("s");
    k++;
    variable.dim[k].name=strdup("eta");
    k++;
    variable.dim[k].name=strdup("xi");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+3)
      );
    k=0;
    sprintf(tmpname,"s_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"eta_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"xi_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  


  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("axis-roms");
  variable.att[k].data=strdup("SETAXI");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("ZTX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
 if(l==0) {
     tmpname =  strdup("level lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(3*l + 22)
      );
    sprintf(tmpname,"levels_%s lat_%s lon_%s",
            grid.name,grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  memcpy(variable.att[k].data,&scale,sizeof(float));
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}

/*----------------------------------------------------------------------------*/

cdfvar_t poc_standardvariable_xietast(char *name, float mask,char *units,float scale,float offset,
                                   char *standardname,char *longname,char *shortname,pocgrd_t grid)

/*----------------------------------------------------------------------------*/
{

  cdfvar_t variable,v1;
  int k,l;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(4);

  variable.initatt(11);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("time");
    k++;
    variable.dim[k].name=strdup("s");
    k++;
    variable.dim[k].name=strdup("eta");
    k++;
    variable.dim[k].name=strdup("xi");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+3)
      );
    k=0;
    
    sprintf(tmpname,"time");
    variable.dim[k].name=strdup(tmpname);
    k++;
                        
    sprintf(tmpname,"s_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"eta_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"xi_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  
  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("axis-roms");
  variable.att[k].data=strdup("TSETAX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TZYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("time levels lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(3*l + 23)
      );
    sprintf(tmpname,"time levels_%s lat_%s lon_%s",
            grid.name,grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}


/*----------------------------------------------------------------------------*/

cdfvar_t poc_standardvariable_xietawt(char *name, float mask,char *units,float scale,float offset,
                                   char *standardname,char *longname,char *shortname,pocgrd_t grid)

/*----------------------------------------------------------------------------*/
{

  cdfvar_t variable,v1;
  int k,l;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(4);

  variable.initatt(11);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("time");
    k++;
    variable.dim[k].name=strdup("w");
    k++;
    variable.dim[k].name=strdup("eta");
    k++;
    variable.dim[k].name=strdup("xi");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+3)
      );
    k=0;
    
    sprintf(tmpname,"time");
    variable.dim[k].name=strdup(tmpname);
    k++;
                        
    sprintf(tmpname,"w");
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"eta_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"xi_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  
  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("axis-roms");
  variable.att[k].data=strdup("TWETAXI");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TWYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("time freq lat lon");
    }
  else {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 19)
      );
    sprintf(tmpname,"time freq lat_%s lon_%s",
            grid.name,grid.name,grid.name);
    }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}
/*----------------------------------------------------------------------------*/

cdfvar_t poc_standardvariable_xietasw(char *name, float mask,char *units,float scale,float offset,
                                   char *standardname,char *longname,char *shortname,pocgrd_t grid)

/*----------------------------------------------------------------------------*/
{

  cdfvar_t variable,v1;
  int k,l;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(4);

  variable.initatt(11);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("w");
    k++;
    variable.dim[k].name=strdup("s");
    k++;
    variable.dim[k].name=strdup("eta");
    k++;
    variable.dim[k].name=strdup("x");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+3)
      );
    k=0;
    sprintf(tmpname,"w");
    variable.dim[k].name=strdup(tmpname);

    k++;
    sprintf(tmpname,"s_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
 
    k++;
    sprintf(tmpname,"eta_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"xi_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  
  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("axis-roms");
  variable.att[k].data=strdup("TWETAXI");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TWYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("time level lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(2*l + 20)
      );
    sprintf(tmpname,"time level lat_%s lon_%s",
            grid.name,grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}

/*----------------------------------------------------------------------------*/

cdfvar_t poc_standardvariable_xietaswt(char *name, float mask,char *units,float scale,float offset,
                                   char *standardname,char *longname,char *shortname,pocgrd_t grid)

/*----------------------------------------------------------------------------*/
{
  cdfvar_t variable,v1;
  int k,l;
  char *tmpname=NULL;

  variable.id=-1;
  variable.type=NC_FLOAT;

  variable.initdim(5);

  variable.initatt(11);

  variable.name=strdup(name);

  l = strlen(grid.name);
  if(l==0) {
    k=0;
    variable.dim[k].name=strdup("time");
    k++;
    variable.dim[k].name=strdup("w");
    k++;
    variable.dim[k].name=strdup("s");
    k++;
    variable.dim[k].name=strdup("eta");
    k++;
    variable.dim[k].name=strdup("xi");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(l+3)
      );
    k=0;
    
    sprintf(tmpname,"time");
    variable.dim[k].name=strdup(tmpname);
    k++;
                        
    sprintf(tmpname,"w");
    variable.dim[k].name=strdup(tmpname);
    k++;

    sprintf(tmpname,"s_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    k++;
    sprintf(tmpname,"eta_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    
    k++;
    sprintf(tmpname,"xi_%s",grid.name);
    variable.dim[k].name=strdup(tmpname);
    free(tmpname);
  }
  
  k=0;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("units");
  variable.att[k].data=strdup(units);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("long_name");
  variable.att[k].data=strdup(longname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("standard_name");
  variable.att[k].data=strdup(standardname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("short_name");
  variable.att[k].data=strdup(shortname);
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("axis-roms");
  variable.att[k].data=strdup("TWSETAXI");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
/* *------------------------------------------------------------------------
  attribute name changed to comply with CF standard*/
//  variable.att[k].name=strdup("axis");
  variable.att[k].name=strdup("content");
  variable.att[k].data=strdup("TWZYX");
  variable.att[k].length=strlen( variable.att[k].data);

  k++;
  variable.att[k].type=NC_CHAR;
  variable.att[k].name=strdup("associate");
  if(l==0) {
     tmpname =  strdup("time freq levels lat lon");
  }
  else {
    exitIfNull(
      tmpname=(char *) malloc(3*l + 27)
      );
    sprintf(tmpname,"time freq levels_%s lat_%s lon_%s",
            grid.name,grid.name,grid.name);
  }
  variable.att[k].data=strdup(tmpname);
  variable.att[k].length=strlen( variable.att[k].data);
  free(tmpname);

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("missing_value");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("_FillValue");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &mask +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("scale_factor");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &scale +l);
  variable.att[k].length=1;

  k++;
  variable.att[k].type=NC_FLOAT;
  variable.att[k].name=strdup("add_offset");
  exitIfNull(
    variable.att[k].data=(char *) malloc(sizeof(float))
    );
  for(l=0;l<sizeof(float);l++) variable.att[k].data[l]=(char) *((char *) &offset +l);
  variable.att[k].length=1;

   return(variable);
}
