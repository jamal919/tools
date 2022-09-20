
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief sealevel time serie input/output
*/
/*----------------------------------------------------------------------------*/

#include "config.h"

#include <stdio.h>
#include <string.h>
#include <fstream>
#include <iostream>
#include <string>
#include <map>

using namespace std;

#include "tools-structures.h"

#include "poc-netcdf-data.hpp"
#include "poc-time.h"
#include "mgr-converter.h"
#include "statistic.h"
#include "functions.h"
#include "spectrum.h"
#include "filter.h"
#include "geo.h"

using namespace std;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void line_err_msg(const string & line,int cnt,int recordI)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string msg=replace(line,"\r","",-1);
  STDOUT_BASE_LINE/*avoid perlReplace.pl*/("#ERROR %d fields while reading record %d:"+msg+"\n",cnt,recordI+1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void scan_hms_z_error_trap(const string & line,int *cnt,double *second,double *z,int recordI)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if( *cnt == 6 ) {
    *z=*second;
    *second=0.0;
    (*cnt)++;
    }
  
  if( *cnt != 7 && recordI>=0) {
    line_err_msg(line,*cnt,recordI);
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int scan_ymd_hms_z(const string & line,int *year,int *month,int *day,int *hour,int *minute,double *second,double *z,int recordI=-1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  scan the most generic and standard ascii line of tide-gauge time series
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int cnt;
  
  cnt = sscanf(line.c_str(), "%4d%*[-/ ]%2d%*[-/ ]%2d %2d%*[: ]%2d%*[: ]%lf %lf", year, month, day, hour, minute, second, z);
  scan_hms_z_error_trap(line,&cnt,second,z,recordI);
  
  return cnt;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int scan_dmy_hms_z(const string & line,int *year,int *month,int *day,int *hour,int *minute,double *second,double *z,int recordI=-1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  scan the most generic and standard ascii line of tide-gauge time series
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int cnt;
  
  cnt = sscanf(line.c_str(), "%2d%*[-/ ]%2d%*[-/ ]%4d %2d%*[: ]%2d%*[: ]%lf %lf", day, month, year, hour, minute, second, z);
  scan_hms_z_error_trap(line,&cnt,second,z,recordI);
  
  return cnt;
}




/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadBODC(const char *filename, tseries_t **serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

BODC Request Format Std. V1.0           Headers=  23 Data Cycles=  70128 BODC QC
Series:  530074 Inv:SFPG 10457           Produced:2002/08/21
Id:DPC9496A.I4D       United Kingdom    Proudman Oceanographic Lab., Bidston, UK
 58d21.8mS056d21.3mW                     Start:19941119005157 End:19961118123646
Depth: floor 3776.0 sensor   -9.0        Nom. sample int.:  15.00     minutes
    6 Parameters included:
Parameter  f  P  Q     Absent Data Value  Minimum Value  Maximum Value     Units
PPSSPS01   Y 30 37              -999.000         -1.061          0.893  decibars
Relative pressure (Bottom mounted pressure sensor)
PPSSPS02   Y 40 47              -999.000         -1.225          0.956  decibars
Relative pressure (Bottom mounted pressure sensor)
PPSSPS03   Y 50 57              -999.000         -1.135          0.902  decibars
Relative pressure (Bottom mounted pressure sensor)
TEMPPR01   Y 60 67                -9.000          0.180          0.580     deg C
Sea temperature (Unspecified temperature probe)
TEMPPR02   Y 70 77                -9.000         -0.320          0.060     deg C
Sea temperature (Unspecified temperature probe)
TEMPPR03   Y 80 87                -9.000          0.550          0.940     deg C
Sea temperature (Unspecified temperature probe)
    1 FORTRAN format record:
(I7,A20,A1,1X,F8.3,A1,1X,F8.3,A1,1X,F8.3,A1,1X,F8.3,A1,1X,F8.3,A1,1X,F8.3,A1)
  Cycle    Date      Time    PPSSPS01  PPSSPS02  PPSSPS03  TEMPPR01  TEMPPR02  TEMPPR03
 Number yyyy mm dd hh mi ssf         f         f         f         f         f         f
      1 1994/11/19 00.51.58     0.210     0.327     0.148     0.580    -0.050     0.830
      2 1994/11/19 01.06.58     0.176     0.277     0.106     0.450    -0.130     0.760
      3 1994/11/19 01.21.58     0.135     0.241     0.072     0.400    -0.150     0.740
-----------------------------------------------------------------------------*/
{
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));
  
  int l, m, n, dum, cnt = 0, nitems;
  int nparameters, nheaders,P,Q,ntargets,*target;
  std::string line, sub;
  double h,*z,*zmin,*zmax,*mask;
  int   lon, lat;
  double lonmin,latmin,factor;
  char c1,c2;
  char s1[128],s2[128], **units, **names;
  size_t pos;
  
  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  int second = 0;
  
  std::getline( input, line );

  pos=line.find ("Headers=", 0);
  sub = line.substr (pos);
  nitems=sscanf(sub.c_str(),"%s %d", s1,&nheaders);
  
  for(k=0;k<2;k++) std::getline( input, line );

  std::getline( input, line );
  nitems=sscanf(line.c_str(), "%2dd%4lfm%c%3dd%4lfm%c", &lat,&latmin,&c1,&lon,&lonmin,&c2);
  switch (c1) {
    case 'N':
      mooring->lat=lat+latmin/60.0;
      break;
    case 'S':
      mooring->lat=-(lat+latmin/60.0);
      break;
    default:
      TRAP_ERR_EXIT(-1, "format error");
    }
  switch (c2) {
    case 'E':
      mooring->lon=lon+lonmin/60.0;
      break;
    case 'W':
      mooring->lon=-(lon+lonmin/60.0);
      break;
    default:
      TRAP_ERR_EXIT(-1, "format error");
    }
  std::getline( input, line );
  nitems=sscanf(line.c_str(), "%s %s %lf", s1, s2, &h);
  mooring->depth=-h;

  std::getline( input, line );
  nitems=sscanf(line.c_str(), "%d", &nparameters);
  
  z=   new double[nparameters];
  zmin=new double[nparameters];
  zmax=new double[nparameters];
  mask=new double[nparameters];
  
  names=new char *[nparameters];
  units=new char *[nparameters];
  
  for(k=0;k<nparameters;k++) {
    names[k]=new char [64];
    units[k]=new char [64];
    }
  
  std::getline( input, line );
  for(k=0;k<nparameters;k++) {
    std::getline( input, line );
    cnt = sscanf(line.c_str(), "%s %s %d %d %lf %lf %lf %s", names[k], s1, &P, &Q, &mask[k], &z[k], &z[k], units[k]);
    std::getline( input, line );
    }
  
  for(k=0;k<4;k++) std::getline( input, line );

  int nrecords = mgr_line_count(filename)-nheaders;
  double *t    = new double[nrecords];
  
  ntargets=0;
  target=new int[nparameters];
  for(k=0;k<nparameters;k++) {
    if(names[k][0]=='P') {
      target[ntargets]=k;
      ntargets++;
      }
    }
  
  *serie=new tseries_t[ntargets];
  for(k=0;k<ntargets;k++) {
    (*serie)[k].x=new double*[1];
    (*serie)[k].x[0]=new double[nrecords];
    (*serie)[k].t=t;
    (*serie)[k].n=nrecords;
    (*serie)[k].nparam=1;
    (*serie)[k].mask=mask[target[k]];
    }
  
  if(strcmp(units[target[0]],"mbars")==0) factor=10.;
  else  factor=1000.;
  n=0;
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
//    cnt = sscanf(line.c_str(), "%7d %4d/%2d/%2d %2d.%2d.%2d%c %8lf%c", &dum, &year, &month, &day, &hour, &minute, &second, &c1, &z[0], &c2);
    std::istringstream data(line);
//    data >> dum >> year >> c1 >> month >> c1 >> day >> hour >> c1 >> minute >> c1 >> second >> c1;
    data >> dum >> year >> c1 >> month >> c1 >> day >> hour >> c1 >> minute >> c1 >> second;
    for(k=0;k<nparameters;k++) {
//      data >> z[k] >> c1;
      data >> s1;
      l=strlen(s1);
      if(s1[l]=='N') s1[l]=0;
      if(s1[l]=='f') s1[l]=0;
      sscanf(s1, "%lf", &z[k]);
      }
//     if( cnt != 10) {
//       line_err_msg(line,cnt,i);
//       continue;
//       }
//     if( c2 != ' ') {
//       cout << "#FLAG: reading  "<< endl << line << endl;
//       continue;
//       }
    for(k=0;k<ntargets;k++) {
      m=target[k];
      if(z[m]!=mask[m]) {
        (*serie)[k].x[0][n]=factor*z[m]/1025.0;
        }
      else {
        (*serie)[k].x[0][n]=mask[m];
        }
      }
    t[n] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    n++;
    }

  for(k=0;k<ntargets;k++) (*serie)[k].first=poctime_getdatecnes(t[0],'d');

  return(ntargets);
  
}


int mgr_loadPuertos2(const char *filename, tseries_t *serie);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadPuertos(const char *filename, tseries_t **serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

 RED DE MAREOGRAFOS : REDMAR

 CODIGO DE LA ESTACION : 3109

 Los NIVELES DEL MAR HORARIOS (Niv_H) se obtienen tras aplicar un
 filtro digital centrado a los datos observados cada 5 minutos.
 La aplicacion de dicho filtro elimina cualquier componente
 de energia con periodo inferior a 1 hora.

 CANALES DE PROCEDENCIA

 Los datos de este mareografo esta compuestos por dos
 subconjuntos de datos diferenciados por el canal de
 llegada de la informacion y por el tipo de tratamiento
 y control de calidad. Estos subconjunto son

  - A)  Datos HISTORICOS
  - B)  Datos recibidos en TIEMPO REAL

 El conjunto A) esta formado por datos que han pasado
 un control de calidad exhaustivo que implica eliminacion
 de valores anomales, asi como control de la estabilidad
 de las referencias o de los desfases de reloj

 El conjunto B) esta formado por datos recibidos en
 TIEMPO REAL sobre los que se ha realizado un control
 de calidad que elimina valores fuera de rango, saltos y
 estabilizaciones anomalas. No se incluyen correciones
 debidas a derivas de reloj o de referencia.

 Para el conjunto A) los valores de marea astronomica
 de cada anyo se calculan con las constantes armonicas
 propias de cada anyo. Sin embargo los niveles de marea
 astronomica de los datos del conjunto B) estan calculados
 a traves de  constantes armonicas promediadas sobre todos
 los anyos disponibles
 LISTADO DE PARAMETROS

 Fecha  GMT  con formato anyo, mes, dia, hora
     AA     : Anyo
     MM     : Mes
     DD     : Dia
     HH     : Hora

 Datos de NIVEL DEL MAR HORARIOS (serie filtrada )

 Niv_H  : Nivel del Mar                                      (cm)
 Mar_H  : Marea   o Componente Astronomica                   (cm)
 Res_H  : Residuo o Componente Meteorologica (Niv_H - Mar_H) (cm)

 Datos de NIVEL DEL MAR CADA 5 MINUTOS (serie observada)

 Niv_00 : Nivel del Mar  a los  00 minutos (cm)
 Niv_05 :  ''       ''          05   ''     ''
 Niv_10 :  ''       ''          10   ''     ''
 Niv_15 :  ''       ''          15   ''     ''
 Niv_20 :  ''       ''          20   ''     ''
 Niv_25 :  ''       ''          25   ''     ''
 Niv_30 :  ''       ''          30   ''     ''
 Niv_35 :  ''       ''          35   ''     ''
 Niv_40 :  ''       ''          40   ''     ''
 Niv_45 :  ''       ''          45   ''     ''
 Niv_50 :  ''       ''          50   ''     ''
 Niv_55 :  ''       ''          55   ''     ''

 Pro    : Especifica el conjunto de datos al que pertenece
           el registro. Puede tomar los siguiente valores.
   1  -  A)  Datos HISTORICOS
   3  -  B)  Datos de TIEMPO REAL


 DATO NULO  : Es representado por -9999

 LISTADO DE DATOS

YY   MM DD HH  Niv_H  Mar_H  Res_H Niv_00 Niv_05 Niv_10 Niv_15 Niv_20 Niv_25 Niv_30 Niv_35 Niv_40 Niv_45 Niv_50 Niv_55 Pro

2008 02 01 00  -9999    329  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999   1
2008 02 01 01  -9999    302  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999   1
2008 02 01 02  -9999    270  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999  -9999   1
-----------------------------------------------------------------------------*/
{
  cout << endl << "-------- starting conversion --------" << endl << endl;

  int m, n;
  int nparameters, nheaders,ntargets,*target;
  std::string line, sub;
  double *z,*zmin,*zmax,*mask;
  double factor;
  char **units, **names;
  size_t pos;
  
  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0/*, minute = 0*/;
//  int second = 0;
  
  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));

  mooring->lon=0;
  mooring->lat=0;
  mooring->depth=0;
  mooring->name=strdup("undocumented");
  
  nheaders=0;
  pos=string::npos;
  while (pos==string::npos) {
    std::getline( input, line );
    pos=line.find ("LISTADO DE DATOS", 0);
    nheaders++;
    }
  
  std::getline( input, line );nheaders++;
  std::getline( input, line );nheaders++;
  pos=line.find (" mm ", 0);
  if(pos!=string::npos){
    *serie=new tseries_t[1];
    mgr_loadPuertos2(filename, *serie);
    return 1;
    }
  std::getline( input, line );nheaders++;
  
  nparameters=15;
  z=   new double[nparameters];
  zmin=new double[nparameters];
  zmax=new double[nparameters];
  mask=new double[nparameters];
  
  names=new char *[nparameters];
  units=new char *[nparameters];
  
  for(k=0;k<nparameters;k++) {
    names[k]=new char [64];
    units[k]=new char [64];
    }
  
  int nrecords = mgr_line_count(filename)-nheaders;
  double *t    = new double[nrecords];
  
  ntargets=1;
  target=new int[nparameters];
  
  target[0]=1;
  
  *serie=new tseries_t[ntargets];
  for(k=0;k<ntargets;k++) {
    (*serie)[k].x=new double*[1];
    (*serie)[k].x[0]=new double[nrecords];
    (*serie)[k].t=t;
    (*serie)[k].n=nrecords;
    (*serie)[k].nparam=1;
    (*serie)[k].mask=-9999;
    }
  
  factor=1.e-02;
  
  n=0;
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    std::istringstream data(line);
    data >> year >>  month  >> day >> hour;
    for(k=0;k<nparameters;k++) {
      data >> z[k];
      }
    for(k=0;k<ntargets;k++) {
      m=target[k];
      if(z[m]!=-9999) {
        (*serie)[k].x[0][n]=factor*z[m];
        }
      else {
        (*serie)[k].x[0][n]=mask[m];
        }
      }
    t[n] = julian_day(month, day, year)  - CNES0jd + hour /24.;
//     t[n] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    n++;
    }

  for(k=0;k<ntargets;k++) (*serie)[k].first=poctime_getdatecnes(t[0],'d');

  return(ntargets);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadPuertos2(const char *filename, tseries_t *serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

 RED DE MAREOGRAFOS : REDMAR

 CODIGO DE LA ESTACION : 3219

 CANALES DE PROCEDENCIA

 Los datos de este fichero esta compuestos por  datos
 HISTORICOS. Este conjunto esta formado por datos que
 han pasado un control de calidad exhaustivo que implica
 la eliminacion de valores anomales, asi como control
 la estabilidad de las referencias o de los desfases
 de reloj

 LISTADO DE PARAMETROS

 Fecha  GMT  con formato anyo, mes, dia, hora
     AA     : Anyo
     MM     : Mes
     DD     : Dia
     HH     : Hora
     MM     : Minuto

    Nivel  : Nivel del Mar   (cm)


 DATO NULO  : Es representado por -9999

 LISTADO DE DATOS

YY   MM DD HH  mm  Nivel

1992 07 01 00 00     224
1992 07 01 00 05     232
1992 07 01 00 10     240
-----------------------------------------------------------------------------*/
{
  std::string line;
//   cout << endl << "-------- starting conversion --------" << endl << endl;
  
  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));
  
/*----------------------------------------------------------------------------
  scan header */
  int nheaders;
  size_t pos;
  
  nheaders=0;
  pos=string::npos;
  while (pos==string::npos) {
    std::getline( input, line );
    nheaders++;
    pos=line.find ("LISTADO DE DATOS", 0);
    }
  
  std::getline( input, line );nheaders++;
  std::getline( input, line );nheaders++;
  pos=line.find (" mm ", 0);
  if(pos==string::npos)
    TRAP_ERR_EXIT(ENOEXEC,"programming error : %s is not a PUERTOS version 2 format file\n",filename);
  
  std::getline( input, line );nheaders++;
  
/*----------------------------------------------------------------------------
  read data */
  int nrecords = mgr_line_count(filename)-nheaders;
  
  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];
  const double mask=-9999;
  
  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;
  
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    int cnt = scan_ymd_hms_z(line, &year, &month, &day, &hour, &minute, &second, &z, i);
    if( cnt != 7) continue;
    if(z!=mask)
      z*=0.01;
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    k++;
    }
  
  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->n=k;
  serie->nparam=1;
  serie->mask=mask;
  serie->first=poctime_getdatecnes(t[0],'d');
  
  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadDART(const char *filename, tseries_t **serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

48.305 N 174.212 W
#YY  MM DD hh mm ss T   HEIGHT
#yr  mo dy hr mn  s -        m
2011 01 01 00 00 00 1 5450.718
2011 01 01 00 15 00 1 5450.646
2011 01 01 00 30 00 1 5450.577
2011 01 01 00 45 00 1 5450.505

-----------------------------------------------------------------------------*/
{
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));
  
  int m, n, nitems;
  int flag, nparameters, nheaders,ntargets,*target;
  std::string line, sub;
  double *z,*zmin,*zmax,*mask;
  double factor;
  char c1,c2;
  char **units, **names;
  
  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  int second = 0;
  
  std::getline( input, line );
  nitems=sscanf(line.c_str(), "%lf %c %lf %c", &mooring->lat,&c1,&mooring->lon,&c2);
  switch (c1) {
    case 'N':
      break;
    case 'S':
      mooring->lat*=-1;
      break;
    default:
      TRAP_ERR_EXIT(-1, "format error");
    }
  switch (c2) {
    case 'E':
      break;
    case 'W':
      mooring->lon*=-1;
      break;
    default:
      TRAP_ERR_EXIT(-1, "format error");
    }
  std::getline( input, line );
  std::getline( input, line );
  
  nparameters=1;
  
  z=   new double[nparameters];
  zmin=new double[nparameters];
  zmax=new double[nparameters];
  mask=new double[nparameters];
  
  names=new char *[nparameters];
  units=new char *[nparameters];
  
  for(k=0;k<nparameters;k++) {
    names[k]=new char [64];
    units[k]=new char [64];
    }
  nheaders=3;
  int nrecords = mgr_line_count(filename)-nheaders;
  double *t    = new double[nrecords];
  
  ntargets=0;
  target=new int[nparameters];
  for(k=0;k<nparameters;k++) {
    target[ntargets]=k;
    ntargets++;
    }
  
  *serie=new tseries_t[ntargets];
  for(k=0;k<ntargets;k++) {
    (*serie)[k].x=new double*[1];
    (*serie)[k].x[0]=new double[nrecords];
    (*serie)[k].t=t;
    (*serie)[k].n=nrecords;
    (*serie)[k].nparam=1;
    (*serie)[k].mask=mask[target[k]];
    }
  
  factor=1.;
  
  n=0;
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    std::istringstream data(line);
    data >> year >> month >> day >> hour >> minute >> second;
    for(k=0;k<nparameters;k++) {
      data >> flag >> z[k];
      }
    for(k=0;k<ntargets;k++) {
      m=target[k];
      if(z[m]!=mask[m]) {
        (*serie)[k].x[0][n]=factor*z[m];
        }
      else {
        (*serie)[k].x[0][n]=mask[m];
        }
      }
    t[n] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    n++;
    }

  for(k=0;k<ntargets;k++) (*serie)[k].first=poctime_getdatecnes(t[0],'d');
  
  mooring->depth=(*serie)[0].x[0][0];
  
  return(ntargets);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadRMN(char *filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  Load RMN data

  1992-01-01 00:00:00.0,-0.28

-----------------------------------------------------------------------------*/
{
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));


  /* skip header */
  int cnt = 0;
  std::string line;
  do{
    std::getline( input, line );
    cnt++;
    }
  while( line.compare(0, 12, "Data/Ora,LIV") != 0);

  /* read data */
  int nlines = mgr_line_count(filename);

  double *t        = new double[nlines - cnt];
  double *sealevel = new double[nlines - cnt];

  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;

  for(int i = 0; i < nlines - cnt; i++){
    std::getline( input, line );
    cnt = line.find_first_of(',');
    line.insert(cnt, " ");
    line.insert(cnt + 2, " ");
    if( ( cnt = sscanf(line.c_str(), "%4d-%02d-%02d %02d:%02d:%lf %lf",
                       &year, &month, &day, &hour, &minute, &second, &(sealevel[k])) ) != 7) {
      cout << "#ERROR: read "<< endl << line << endl;
      continue;
      }
    t[k++] = ( julian_day(month, day, year) - CNES0jd )
           + hour /24.
           + minute / (24. * 60.0);
    }

  delete [] t;
  delete [] sealevel;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadRadar(const char *filename, tseries_t **serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  Load RADAR-HF netcdf data

-----------------------------------------------------------------------------*/
{
  int n,ncid,keep,status,verbose=0;
  int lon_id,lat_id,u_id,v_id,t_id,mask_id;
  int nx,ny,nt,np;
  double *u, *v, *mask,spec=9999.;
  double *t,*lon,*lat;
  cdfvar_t info;
  cdfgbl_t global;
  
  cout << endl << "-------- starting conversion --------" << endl << endl;
  status= cdf_globalinfo(filename,&global,verbose);

//  status=cdf_varinfo(filename,varid,&info);
  
  lon_id= cdf_identify(global, "Lon");
  lat_id= cdf_identify(global, "Lat");
  t_id= cdf_identify(global, "T");
  u_id= cdf_identify(global, "U");
  v_id= cdf_identify(global, "V");
  mask_id= cdf_identify(global, "mask");
  
  status=cdf_varinfo(filename,lon_id,&info);
  nx=info.dim[0].length;

  status=cdf_varinfo(filename,lat_id,&info);
  ny=info.dim[0].length;

  status=cdf_varinfo(filename,t_id,&info);
  nt=info.dim[0].length;

  lon=new double[nx];
  lat=new double[ny];
  t=new double[nt];
  
  status=nc_open(filename,NC_NOWRITE,&ncid);
  status=nc_get_var_double (ncid,lon_id,lon);
  status=nc_get_var_double (ncid,lat_id,lat);
  status=nc_get_var_double (ncid,t_id,t);
  
  u=new double[nx*ny*nt];
  v=new double[nx*ny*nt];
  mask=new double[nx*ny];
  status=nc_get_var_double (ncid,u_id,u);
  status=nc_get_var_double (ncid,v_id,v);
  status=nc_get_var_double (ncid,mask_id,mask);
  
  np=0;
  for(n=0;n<nx*ny;n++) {
    if(mask[n]!=0) np++;
    }
  
  for(n=0;n<nt*nx*ny;n++) {
    if(isnan(u[n])!=0) u[n]=spec;
    if(isnan(v[n])!=0) v[n]=spec;
    }
  
  for(n=0;n<nx*ny;n++) {
    keep=0;
    for(size_t k=0;k<nt;k++) {
      if(u[k*nx*ny+n]!=spec) {
        keep=1;
        break;
        }
      if(keep==0) mask[n]=0;
      }
    }

  np=0;
  for(n=0;n<nx*ny;n++) {
    if(mask[n]!=0) np++;
    }
  
//   double tt=mjd(1,1,2008)+211;
//   tt=julian_day(1,1,2008)+212-julian_day(1,1,1)-365;
//
//   *serie=new tseries_t[np];

//   for(n=0;n<nt;n++) {
//     t[n]-=julian_day(1,1,1);
// //    t[n]-=mjd(1,1,1950);
//     }

  *serie=new tseries_t[np];

  np=0;
  for(n=0;n<nx*ny;n++) {
    if(mask[n]!=0){
      (*serie)[np].x=new double*[2];
      (*serie)[np].x[0]=new double[nt];
      (*serie)[np].x[1]=new double[nt];
      for(size_t k=0;k<nt;k++) {
        (*serie)[np].x[0][k]=u[k*nx*ny+n];
        (*serie)[np].x[1][k]=v[k*nx*ny+n];
        }
      (*serie)[np].n=nt;
      (*serie)[np].t=t;
/* *----------------------------------------------------------------------------
      m=i*ny+j, j=m%nx, i=m/ny*/
      (*serie)[np].lon=lon[n/ny];
      (*serie)[np].lat=lat[n%ny];
      (*serie)[np].mask=spec;
      (*serie)[np].nparam=2;
      np++;
      }
    }

  return(np);
}


#define GLOBAL  0
#define YEARLY  1
#define MONTHLY 2
#define WEEKLY  3
#define DAILY   4
#define HOURLY  5
#define CUSTOM  6

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_name(date_t actual, const char *name_template, const char *varname, char **filename, int *packing)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  status;
  char dummy[256], *pointer;
  date_t cdf_reference;
  FILE *out;

/* *----------------------------------------------------------------------------
  Reconstruct the input file name*/
  (*filename) = new char[strlen(name_template) + 256];
  sprintf((*filename), "%s", name_template);

  out = fopen(*filename, "r");
  status = (out == NULL);

  switch (status) {
    case 0:
/* *----------------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      *packing=GLOBAL;
      break;

    default:
/*-----------------------------------------------------------------------------
      use format convention information*/
/* *----------------------------------------------------------------------------
      substitute YYYY with current year*/
      pointer = strstr((*filename), "YYYY");
      if(pointer != NULL) {
        sprintf(dummy, "%4d", actual.year);
        strncpy(pointer, dummy, 4);
        *packing=YEARLY;
        }
/* *----------------------------------------------------------------------------
      substitute MM with current month*/
      pointer = strstr((*filename), "MM");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.month);
        strncpy(pointer, dummy, 2);
        *packing=MONTHLY;
        }
/* *----------------------------------------------------------------------------
      substitute DD with current day*/
      pointer = strstr((*filename), "DD");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.day);
        strncpy(pointer, dummy, 2);
        *packing=DAILY;
        }
/* *----------------------------------------------------------------------------
      substitute DD with current day*/
      pointer = strstr((*filename), "NNN");
      if(pointer != NULL) {
        int day=day_in_year(actual);
        sprintf(dummy, "%3.3d", day);
        strncpy(pointer, dummy, 3);
        *packing=DAILY;
        }
/* *----------------------------------------------------------------------------
      substitute HH with current hour*/
      pointer = strstr((*filename), "HH");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", (int) floor(actual.second/3600.));
        strncpy(pointer, dummy, 2);
        *packing=HOURLY;
        }
/* *----------------------------------------------------------------------------
      substitute MN with current hour*/
      pointer = strstr((*filename), "MN");
      if(pointer != NULL) {
        double hour=floor(actual.second/3600.);
        sprintf(dummy, "%2.2d", (int) floor((actual.second-3600.*hour)/60.));
        strncpy(pointer, dummy, 2);
        *packing=HOURLY;
        }
// /* *----------------------------------------------------------------------------
//       lookup file fitting ???? with minute and seconds index*/
//       pointer = strstr((*filename), "????");
//       if(pointer != NULL) {
//         int i,j;
//         for(j=0;j<60;j++) {
//           for(i=0;i<60;i++) {
//             sprintf(dummy, "%2.2d%2.2d", j,i);
//             strncpy(pointer, dummy, 4);
//             status= get_file_size(*filename, 0);
// /* *----------------------------------------------------------------------------
//             status=0 if file found*/
//             if(status==0) break;
//             }
//           if(status==0) break;
//           }
//         }
/* *----------------------------------------------------------------------------
      substitute VARNAME*/
      pointer = strstr((*filename), "VARNAME");
      if(pointer != NULL) {
        sprintf(dummy, "%s", varname);
        strncpy(pointer, dummy, strlen(varname));
        }
      break;

      break;
    }

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadRadarRaw(const char *filename, tseries_t **serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Load RADAR-HF netcdf data

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k, l, i, n, m, nitems, keep, status;
  int nx,ny,nt,np;
  double *u, *v, *mask,spec=9999.;
  double *time, *lon, *lat,t0, *azimuth, *radius;
  date_t start[2], end;
  FILE *in;
  
  cout << endl << "-------- starting conversion --------" << endl << endl;
  
  ny=71;
  nx=100;
  nt=70*24*3;

  lon=new double[nx*ny];
  lat=new double[nx*ny];
  time=new double[nt];
    
  u=new double[nx*ny*nt];
  v=new double[nx*ny*nt];

  for(n=0;n<nt*nx*ny;n++) {
    u[n]=spec;
    v[n]=spec;
    }
  
  mask=new double[nx*ny];
  
  azimuth=new double[ny];
  radius=new double[nx];
 
  double azimuth_shift[2]={-21.,+20.};
  double ref_lat[2],ref_lon[2];
  double x0,y0;
  
  ref_lon[0]= -4.6667;
  ref_lat[0]= 48.069;
  
  ref_lon[1]= -4.7757;
  ref_lat[1]= 48.5030;

  start[0]=date_t(2007, 8,23,600.);
  start[1]=date_t(2007, 8,23,0.);
  
  end  =date_t(2007,10,27,0.);
  
  t0=cnes_time(start[0],'s');
  
  for(i=0;i<nt;i++) {
    time[i]=t0+i*1200.;
    date_t present = poctime_getdatecnes(time[i], 's');
    char *frame, keyword[1024];
    int packing;
    status=decode_name(present, "YYYYNNNHHMN_bre.coc", 0, &frame, &packing);
    in=fopen(frame,"r");
    if(in==0) continue;
//    fgets(line,1024,in);
    printf("reading %s\n",frame);
    fscanf(in,"%s",keyword);
    for(k=0;k<ny;k++) {
      fscanf(in,"%lf",&azimuth[k]);
      }
   
    for(l=0;l<nx;l++) {
      fscanf(in,"%lf",&radius[l]);
      for(k=0;k<ny;k++) {
//        fgets(line,1024,in);
        n=k*ny+l;
        m=i*nx*ny+n;
        fscanf(in,"%s",keyword);
        nitems=sscanf(keyword,"%lf",&u[m]);
        if(nitems!=1) {
          u[m]=spec;
          }
        if(isnan(u[m])!=0) {
          u[m]=spec;
          }
        }
      }
    fclose(in);
    delete[] frame;
    }
  
  projPJ projection;
  projection=assign_StereoOblique(ref_lat[0],ref_lon[0]);
  
  geo_to_projection(projection, ref_lat[0],ref_lon[0],&x0,&y0);
  
  for(k=0;k<ny;k++) {
    for(l=0;l<nx;l++) {
      double alpha=(azimuth[k]+azimuth_shift[0]+90.)*M_PI/180.0;
      double x=x0+radius[l]*cos(alpha)*1000.;
      double y=y0+radius[l]*sin(alpha)*1000.;
      n=k*ny+l;
      double t,p;
      projection_to_geo(projection, &p, &t,x,y);
      lon[n]=t;
      lat[n]=p;
      }
    }
  
    
//  Brezellec: -4.6667  48.069  -21.00 (Master) x=100.4556; y= 63.2617 km (Ref.: lon0=-6; lat0=47.5)
//  Garchine:  -4.7757  48.5030  20.00 (Slave)  x= 92.2432; y=111.5141 en km
//      cartes axe 0° de bas en haut  (n'est plus valable; 90 deg ajoutes)
//         x0=0; y0=0; theta0=0; lon0=-6; lat0=47.5;  % pt de reference pour la grille en km
//         if i_sta==1;
//           theta0=-21+90;
//           [x0 y0]=lonlat2km(lon0,lat0, -4.6667,  48.069);
//           teta=fliplr(teta')';disp('RETOURNEMENT DES AZIMUTHS POUR BREZELLEC (19/10/09)');
//         else
//           theta0=20+90;
//           [x0 y0]=lonlat2km(lon0,lat0, -4.7757,  48.503);
//         end
//         alphal=pi*(theta0+teta+90)/180.;
//         x=x0+dist*cos(alphal');
//         y=y0+dist*sin(alphal');

  
  for(n=0;n<nx*ny;n++) {
    mask[n]=1;
    keep=0;
    for(size_t k=0;k<nt;k++) {
      if(u[k*nx*ny+n]!=spec) {
        keep=1;
        break;
        }
      }
    if(keep==0) mask[n]=0;
    }

  np=0;
  for(n=0;n<nx*ny;n++) {
    if(mask[n]!=0) np++;
    }
  
  *serie=new tseries_t[np];

  np=0;
  for(k=0;k<ny;k++) {
    for(l=0;l<nx;l++) {
      n=k*nx+l;
      if(mask[n]!=0){
        double alpha=(azimuth[k]+azimuth_shift[0]+90.)*M_PI/180.0;
        (*serie)[np].x=new double*[2];
        (*serie)[np].x[0]=new double[nt];
        (*serie)[np].x[1]=new double[nt];
        for(size_t i=0;i<nt;i++) {
          if(u[i*nx*ny+n]!=spec) {
            double ux=u[i*nx*ny+n]*cos(alpha);
            double uy=u[i*nx*ny+n]*sin(alpha);
            (*serie)[np].x[0][i]=ux;
            (*serie)[np].x[1][i]=uy;
            }
          else {
            (*serie)[np].x[0][i]=spec;
            (*serie)[np].x[1][i]=spec;
            }
          }
        (*serie)[np].n=nt;
        (*serie)[np].t=time;
/* *----------------------------------------------------------------------------
        m=i*ny+j, j=m%nx, i=m/ny*/
        (*serie)[np].lon=lon[n];
        (*serie)[np].lat=lat[n];
        (*serie)[np].mask=spec;
        (*serie)[np].nparam=2;
        np++;
        }
      }
    }

  return(np);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadList2(const char *filename, tseries_t **serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/// Load LIST(2) ascii data
/** \return number of times series */
/*
#-- HEADER -------------------------------------------
# Column 1 : date in day referred to
#            CNES date (01-JAN-1950 00:00:00.0)
# Column 2 : sea surface height (in meters)
#-- HEADER END ---------------------------------------
Number of crossover points :
9
#------------------------------------------------------
Pt  : 0
Lon : 359.447179
Lat : 44.859973
Mes : 201047
1609455600.000000 4.060000
1609455900.000000 3.990000

-----------------------------------------------------------------------------*/
{
  int k, n, nitems;
  int target,ntargets, nrecords;
  FILE *in;
  char   line[1024], name[64], keyword[32], c;
    
  cout << endl << "-------- starting conversion --------" << endl << endl;
  
  in = fopen(filename,"r");
  
  if(in == NULL) {
    cout << "ERROR : can not open : "<< filename << endl;
    return 0;
    }

  do  fgets(line, sizeof(line), in);
    while (strncmp(line,"#-- HEADER END ---------------------------------------",20) != 0 );
  
  fgets(line,sizeof(line),in);
  if( (nitems = fscanf(in, "%d\n", &ntargets)) !=1) {
    cout << "ERROR : can not read : "<< line <<endl;
    return 0;
    }
  
  *serie=new tseries_t[ntargets];
  for(k=0;k<ntargets;k++) {
    fgets(line,sizeof(line),in);
    fgets(line,sizeof(line),in);
    nitems=sscanf(line,"%s %c %d %s",  keyword, &c, &target, name);
    if(nitems==3) strcpy(name,"");
    fscanf(in,"%s %c %lf", keyword, &c, &((*serie)[k].lon));
    fscanf(in,"%s %c %lf", keyword, &c, &((*serie)[k].lat));
    fscanf(in,"%s %c %d",  keyword, &c, &nrecords);
    (*serie)[k].nparam=1;
    (*serie)[k].x=new double*[(*serie)[k].nparam];
    (*serie)[k].x[0]=new double[nrecords];
    (*serie)[k].t=new double[nrecords];
    (*serie)[k].n=nrecords;
    (*serie)[k].mask=9999.;
    for(n=0;n<nrecords;n++) {
      fscanf(in,"%lf %lf", &(*serie)[k].t[n], &(*serie)[k].x[0][n]);
//      (*serie)[k].t[n]/=24.*3600;
      }
    fgets(line,sizeof(line),in);
    (*serie)[k].mooring.lon=(*serie)[k].lon;
    (*serie)[k].mooring.lat=(*serie)[k].lat;
    (*serie)[k].mooring.name=strdup(name);
    }
  
  fclose(in);
  
  return(ntargets);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**
\brief loads a binary file containing Profilers time series.

No header is read (for now ....),
First we have the number of bins and the number ot time samples then bin1 bin2 etc.
Everything is a FLOAT.
The read data are stored in a tseries_t object.

@param [in] filename the file to read
@param [out] serie  a tseries_t object containing the read data and possibly some more information (like sampling, mask value, timezone, first and last date)
@return int nrecords : number of records in the time serie
*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadProfilers(const char *filename, tseries_t **serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  Load Profilers data

-----------------------------------------------------------------------------*/
{
  FILE *f_in;
  float nz, nt;
  float *df,*time;
  int DateTime[6];
  int nbin,ntime;
   
  /* lecture du fichier input */
  f_in = fopen(filename, "r");
  if(f_in==0) {
      fprintf(stderr, "Cannot read the file: %s\n",filename);
      return(EXIT_FAILURE);
    }
  fread( (&nz), sizeof(nz),1,f_in);
  fread((&nt), sizeof(nt),1,f_in);
  ntime=int(nt);
  nbin=int(nz);
  printf("the file contains %d bin and %d time samples\n",nbin,ntime);
  //nbin=nb;
  //ntime=nt;
  
  *serie=new tseries_t[nbin];
  for (int j=0;j<nbin;j++) {
    (*serie)[j].x = new double* [1];
    (*serie)[j].x[0]= new double[ntime];
    (*serie)[j].t=new double[ntime];
    (*serie)[j].lon=j;
    (*serie)[j].lat=j;
    (*serie)[j].n=ntime;
    (*serie)[j].nparam=1;
  }
   for (int i=0; i<ntime;i++) {
      time = new float[6];
      fread(time,sizeof(float),6,f_in);
      for (int j=0;j<6;j++) DateTime[j]= int(time[j]);
      delete[] time;
           
      df = new float[nbin];
      fread(df, sizeof(float),nbin,f_in);
      
      for (int j=0;j<nbin;j++) {
        (*serie)[j].x[0][i] = double(df[j]);
        (*serie)[j].t[i]= ( julian_day(DateTime[1],DateTime[2], DateTime[0]) - CNES0jd )
           + DateTime[3]/24.
           + DateTime[4]/ (24. * 60.0)+ DateTime[5]/ (24. * 60.0* 60.0);
      }
      delete[] df;
     }
   
     //delete[] DateTime;
   
  return(nbin);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadPAB(const char *filename, tseries_t *serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  Load PAB (Port autonome de Bordeaux) data

-----------------------------------------------------------------------------*/
{
  std::string line;
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));

/* *----------------------------------------------------------------------------
  read data */
  int nrecords = mgr_line_count(filename);

  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0;

  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    int cnt = sscanf(line.c_str(), "%d %d %d %d %lf", &year, &month, &day, &hour, &z);
    if( cnt != 5) {
      cout << "#ERROR: read "<< endl << line << endl;
      continue;
      }
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24.;
    k++;
    }
  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->n=k;
  serie->first=poctime_getdatecnes(t[0],'d');

  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadMINH(const char *filename, tseries_t *serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------------

  Load MINH (vietnam) data

-----------------------------------------------------------------------------*/
{
  std::string line;
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));

/* *----------------------------------------------------------------------------
  read data */
  int nrecords = mgr_line_count(filename)-2;

  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0;

  std::getline( input, line );
  std::getline( input, line );
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    int cnt = sscanf(line.c_str(), "%lf %d %d %d %d", &z, &hour, &day, &month, &year);
    if( cnt != 5) {
      cout << "#ERROR: read "<< endl << line << endl;
      continue;
      }
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24.;
    k++;
    }
  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->n=k;
  serie->nparam=1;
  serie->mask=-999;
  serie->first=poctime_getdatecnes(t[0],'d');
  
  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**
\brief loads a file containing a SONEL_HR sealevel time serie.

No header is read (for now ....), each line is formated as : YYYY-MM-DD HH:MM:SS    float_value_of_slv   ( ex : 1995-08-10 11:50:00   3.050 )
The read data are stored in a tseries_t object.

@param [in] filename the file to read
@param [out] serie  a tseries_t object containing the read data and possibly some more information (like sampling, mask value, timezone, first and last date)
@return int nrecords : number of records in the time serie
*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadSONEL_HR(const char *filename, tseries_t *serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* NOTE: make sure the record format is NOT simlar to the example below !!! */
/*----------------------------------------------------------------------------

  Load SONEL_HR data
  
# ker_1995.dat
# 69 27.178 E
# 47 40.233 S
# 185 m
21/11/1995 11:00:00  18134.771
21/11/1995 12:00:00  18121.801
21/11/1995 13:00:00  18115.824
21/11/1995 14:00:00  18119.602

-----------------------------------------------------------------------------*/
{
  std::string line;
  char c,c1,c2;
  double lat,lon,latmin,lonmin;
  int nitems;
  
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));

/*----------------------------------------------------------------------------
  read data */
  int nrecords = mgr_line_count(filename);

  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;

  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    if(line[0]=='#') {
      sscanf(line.c_str(), "%c %s", &c, mooring->name);
      std::getline( input, line );
      nitems=sscanf(line.c_str(), "%c %lf %lf %c", &c, &lon,&lonmin, &c2);
      std::getline( input, line );
      nitems=sscanf(line.c_str(), "%c %lf %lf %c", &c, &lat,&latmin, &c1);
      switch (c1) {
        case 'N':
          mooring->lat=lat+latmin/60.0;
          break;
        case 'S':
          mooring->lat=-(lat+latmin/60.0);
          break;
        default:
          TRAP_ERR_EXIT(-1, "format error");
        }
      switch (c2) {
        case 'E':
          mooring->lon=lon+lonmin/60.0;
          break;
        case 'W':
          mooring->lon=-(lon+lonmin/60.0);
          break;
        default:
          TRAP_ERR_EXIT(-1, "format error");
        }
      std::getline( input, line );
      nitems=sscanf(line.c_str(), "%c %lf", &mooring->depth);
      std::getline( input, line );
      nrecords-=4;
      }
    int cnt = scan_dmy_hms_z(line, &year, &month, &day, &hour, &minute, &second, &z, i);
    if( cnt != 7) continue;
    
/**----------------------------------------------------------------------------
    convert mbars to meters */
//    sealevel[k]=z/100./1.025;
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    k++;
    }

  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->mask=999999.0;
  serie->n=k;
  serie->first=poctime_getdatecnes(t[0],'d');

  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadREFMAR(const char *filename, tseries_t *serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------

  Load REFMAR data

# Station : SAINT-NAZAIRE
# Longitude : -2.201550
# Latitude : 47.266862
# Organisme fournisseur de données : SHOM / GPM de Nantes-Saint Nazaire
# Fuseau horaire : UTC
# Unité : m
# Source 4 : Données horaires validées en temps différé
Date Valeur Source
15/05/1957 23:00:00 1.060 4
16/05/1957 00:00:00 2.050 4

-----------------------------------------------------------------------------*/
{
  std::string line, key;
//   char c,dum[1024];
  char *name=new char[1024];
  double lat,lon;
  int nitems;
  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;

  
  cout << endl << "-------- starting conversion --------" << endl << endl;

/*----------------------------------------------------------------------------
  get datafile size */
  int nrecords = mgr_line_count(filename);

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));
  
  std::getline( input, line );
  
  if(line[0]=='#') {
//     sscanf(line.c_str(), "%c %s %c %s", &c, dum, &c, name);
//     mooring->name=strdup(name);
//     std::getline( input, line );
//     nitems=sscanf(line.c_str(), "%c %s %c %lf", &c, dum, &c, &lon);
//     std::getline( input, line );
//     nitems=sscanf(line.c_str(), "%c %s %c %lf", &c, dum, &c, &lat);
    key=line.substr(line.find(":")+2);
    nitems=sscanf(key.c_str(), "%s", name);
    std::getline( input, line );
    key=line.substr(line.find(":")+2);
    nitems=sscanf(key.c_str(), "%lf",&lon);
    std::getline( input, line );
    key=line.substr(line.find(":")+2);
    nitems=sscanf(key.c_str(), "%lf", &lat);
    mooring->name=strdup(name);
    mooring->lat=lat;
    mooring->lon=lon;
    for(k=0;k<5;k++) std::getline( input, line );
    nrecords-=8;
    }
  else {
    }
  
  if(serie==0) {
    input.close();
    return(0);
    }
  
  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  k = 0;
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    int cnt = sscanf(line.c_str(), "%2d/%2d/%4d%*[ ;]%2d:%2d:%2d%*[ ;]%lf", &day, &month, &year, &hour, &minute, &second, &z);
    if( cnt != 7) {
      line_err_msg(line,cnt,i);
      continue;
      }
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    k++;
    }
  input.close();

  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->mask=999999.0;
  serie->n=k;
  serie->first=poctime_getdatecnes(t[0],'d');

  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadBODC2(const char *filename, tseries_t *serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------

  Load BODC harbour data

Port:              P935
Site:              Portrush
Latitude:          55.20678
Longitude:         -6.65683
Start Date:        1995/01/01 00:00:00
End Date:          2014/06/30 23:45:00
Contributor:       National Oceanography Centre, Liverpool
Datum information: UK Admiralty Chart Datum
Parameter code:    ASLVBG02 = Surface elevation (unspecified datum second sensor) of the water body by bubbler tide gauge
  Cycle    Date      Time    ASLVBG02
 Number yyyy mm dd hh mi ssf         f
     1) 1995/01/01 00:00:00   -99.000N
     2) 1995/01/01 00:15:00   -99.000N

-----------------------------------------------------------------------------*/
{
  std::string line, key;
  char *name=new char[1024];
  double lat,lon;
  int nitems;
  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;
  size_t pos;

  
  cout << endl << "-------- starting conversion --------" << endl << endl;

/*----------------------------------------------------------------------------
  get datafile size */
  int nrecords = mgr_line_count(filename);

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));
  
  while(true) {
    std::getline( input, line );
    pos=line.find("Site:              ");
    if(pos!=string::npos) {
      key=line.substr(pos+19);
      nitems=sscanf(key.c_str(), "%s", name);
      }
    pos=line.find("Latitude:          ");
    if(pos!=string::npos) {
      key=line.substr(pos+19);
      nitems=sscanf(key.c_str(), "%lf", &lat);
      }
    pos=line.find("Longitude:         ");
    if(pos!=string::npos) {
      key=line.substr(pos+19);
      nitems=sscanf(key.c_str(), "%lf", &lon);
      }
    pos=line.find("Number");
    if(pos!=string::npos) {
      break;
      }
    nrecords--;
    }
  mooring->name=strdup(name);
  mooring->lat=lat;
  mooring->lon=lon;
    
  if(serie==0) {
    input.close();
    return(0);
    }
  
  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  k = 0;
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    pos=line.find(")");
    if(pos!=string::npos) {
      line=line.substr(pos+2);
      }
    int cnt = scan_ymd_hms_z(line, &year, &month, &day, &hour, &minute, &second, &z, i);
    if( cnt != 7) continue;
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    k++;
    }
  input.close();

  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->mask=-99.0;
  serie->n=k;
  serie->first=poctime_getdatecnes(t[0],'d');

  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadOceanFisheriesCA(const char *filename, tseries_t *serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------

Station_Name,Koksoak River Entrance
Station_Number,4295
Latitude_Decimal_Degrees,58.526139
Longitude_Decimal_Degrees,68.199694
Datum,CD
Time_zone,UTC
SLEV=Observed Water Level
Obs_date,SLEV(metres)
1952/08/23 06:00,10.11,

-----------------------------------------------------------------------------*/
{
  std::string line, key;
  char *name=new char[1024];
  double lat,lon;
  int nitems;
  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;
  size_t pos;

  
  cout << endl << "-------- starting conversion --------" << endl << endl;

/*----------------------------------------------------------------------------
  get datafile size */
  int nrecords = mgr_line_count(filename);

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));
  
  while(true) {
    std::getline( input, line );
    key="Station_Name,";
    pos=line.find(key);
    if(pos!=string::npos) {
      for(size_t k=0; k< line.length();k++) if (line[k]==' ') line[k]='_';
      key=line.substr(pos+key.length());
      nitems=sscanf(key.c_str(), "%s", name);
      }
    key="Latitude_Decimal_Degrees,";
    pos=line.find(key);
    if(pos!=string::npos) {
      key=line.substr(pos+key.length());
      nitems=sscanf(key.c_str(), "%lf", &lat);
      }
    key="Longitude_Decimal_Degrees,";
    pos=line.find(key);
    if(pos!=string::npos) {
      key=line.substr(pos+key.length());
      nitems=sscanf(key.c_str(), "%lf", &lon);
      lon=-lon;
      }
    key="Obs_date,SLEV(metres)";
    pos=line.find(key);
    if(pos!=string::npos) {
      break;
      }
    nrecords--;
    }
  
  mooring->name=strdup(name);
  mooring->lat=lat;
  mooring->lon=lon;
    
  if(serie==0) {
    input.close();
    return(0);
    }
  
  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  k = 0;
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    pos=line.find(")");
    if(pos!=string::npos) {
      line=line.substr(pos+2);
      }
    second=0;
    int cnt = sscanf(line.c_str(), "%4d/%2d/%2d%*[ ,]%2d:%2d%*[ ,]%lf", &year, &month, &day, &hour, &minute, &z);
    if( cnt != 6) {
      STDOUT_BASE_LINE("#ERROR %d!=6 while reading l.%d:"+line+"\n",cnt,i+1);
      continue;
      }
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    k++;
    }
  input.close();

  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->mask=-99.0;
  serie->n=k;
  serie->first=poctime_getdatecnes(t[0],'d');

  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadUser(const char *filename, tseries_t *serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*
##################################################################################################################################################
# Yong-an Tidal Station
# Description: Tip of the breakwater north of Chinese Petroleum Corp. Yong-an LNG Factory LNG Harbor
# Station coordinate: WGS84(120.198,22.819)
# Code station: 1786
# Data source: Central Weather Bureau
# Reference level: TWVD2001 MSL
# Separator: Tabulation \t
# Column 1: year month day hour minute second Zone Asia/Taipei UTC+8
# Column 2: year month day hour minute second UTC
# Column 3: centisecond
# Column 4: Predicted water level from the harmonic components (millimeter) - FillValue: NaN
# Column 5: Predicted water level from the harmonic components (meter) - FillValue: NaN
# Column 6: Observed water level (millimeter) - FillValue: NaN
# Column 7: Observed water level (meter) - FillValue: NaN
##################################################################################################################################################
20110901000000  20110831160000  00      438     0.438   508     0.508
20110901000600  20110831160600  00      426     0.426   488     0.488
20110901001200  20110831161200  00      413     0.413   496     0.496
20110901001800  20110831161800  00      401     0.401   476     0.476
*/
{
  std::string line, key;
  char *name=new char[1024], *dum=new char[1024];
  double lat,lon;
  int nitems;
  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;
  size_t pos;
  
  cout << endl << "-------- starting conversion --------" << endl << endl;

/*----------------------------------------------------------------------------
  get datafile size */
  int nrecords = mgr_line_count(filename);

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));
  
  while(true) {
    std::getline( input, line );
    key=" Tidal Station";
    pos=line.find(key);
    if(pos!=string::npos) {
      string nameStr=replace(line,key);
      replace(&nameStr,"# ");
      replace(&nameStr," ","_");
      strcpy(name,nameStr.c_str());
      }
    key="Station coordinate: WGS84(";
    pos=line.find(key);
    if(pos!=string::npos) {
      key=line.substr(pos+key.length());
      nitems=sscanf(key.c_str(), "%lf,%lf", &lon, &lat);
      }
    key="Column 7";
    pos=line.find(key);
    if(pos!=string::npos) {
      std::getline( input, line );
      break;
      }
    nrecords--;
    }
  
  mooring->name=strdup(name);
  mooring->lat=lat;
  mooring->lon=lon;
    
  if(serie==0) {
    input.close();
    return(0);
    }
  
  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  k = 0;
  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    key="NaN";
    pos=line.find(key);
    if(pos!=string::npos) {
      continue;
      }
    int cnt = sscanf(line.c_str(), "%s%*[ ,]%4d%2d%2d%2d%2d%2d%*[ ,]%s%*[ ,]%s%*[ ,]%s%*[ ,]%s%*[ ,]%lf", dum, &year, &month, &day, &hour, &minute, &second, dum, dum, dum, dum, &z);
    if( cnt != 12) {
      line_err_msg(line,cnt,i);
      continue;
      }
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    k++;
    }
  input.close();

  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->mask=-99.0;
  serie->n=k;
  serie->first=poctime_getdatecnes(t[0],'d');

  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**
\brief loads a file containing a JMAcolon sealevel time serie.

No header is read (for now ....), each line is formated as : YYYY:MM:DD HH:MM  float_value_of_slv   ( ex : 1995:08:10 11:50   3.050 )
The read data are stored in a tseries_t object.

@param [in] filename the file to read
@param [out] serie  a tseries_t object containing the read data and possibly some more information (like sampling, mask value, timezone, first and last date)
@return int nrecords : number of records in the time serie
*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadDMY(const char *filename, tseries_t *serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Load DMY time format data

------------------------------------------------------------------------------*/
{
  std::string line;
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));

/*----------------------------------------------------------------------------
  read data */
  int nrecords = mgr_line_count(filename);

  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;

  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    if(line[0]=='#') continue;
//     int cnt = sscanf(line.c_str(), "%4d:%2d:%2d %2d:%2d %lf", &year, &month, &day, &hour, &minute, &z);
//     if( cnt != 6) {
//       line_err_msg(line,cnt,i);
//       continue;
//       }
    int cnt = scan_dmy_hms_z(line, &year, &month, &day, &hour, &minute, &second, &z, i);
    if( cnt != 7) {
//       line_err_msg(line,cnt,i);
      continue;
      }
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    k++;
    }

  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->n=k;
  serie->first=poctime_getdatecnes(t[0],'d');

  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mgr_loadYMD(const char *filename, tseries_t *serie)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  2010/01/01 00:00:00  1.68
  2010/01/01 00:05:00  1.74

  or

  2006-01-01 00:00:00   1.020
  2006-01-01 00:05:00   0.970
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  std::string line="";
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));

/*----------------------------------------------------------------------------
  read data */
  int nrecords = mgr_line_count(filename);

  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;

  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    if(line.c_str()[0]=='#') continue;
    int cnt = scan_ymd_hms_z(line, &year, &month, &day, &hour, &minute, &second, &z, i);
    z*=0.01;
    if( cnt != 7) continue;
    sealevel[k]=z;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
    k++;
    }

  serie->x=new double*[1];
  serie->x[0]=sealevel;
  serie->t=t;
  serie->n=k;
  serie->mask=99999;
  serie->first=poctime_getdatecnes(t[0],'d');

  return(nrecords);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int mgr_loadSeine(const char *filename, tseries_t *serie, mooring_t *mooring)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------

  Load GIP Seine-AVAL data

#Petit_Couronne 1.01123506 49.38690261
01 01 2007 00:00:00     738
01 01 2007 00:05:00     740
01 01 2007 00:10:00     738

-----------------------------------------------------------------------------*/
{
  char c;
  double lat,lon;
  char name[1024];

  std::string line;
  cout << endl << "-------- starting conversion --------" << endl << endl;

  std::ifstream input(filename);
  if(input.fail())
    TRAP_ERR_RETURN(0,1,"%s: ifstream(\"%s\") error (%d %s)\n",__func__,filename,errno,strerror(errno));
  
  std::getline( input, line );
  if(line[0]=='#') {
    sscanf(line.c_str(), "%c %s %lf %lf",&c,name, &lon,&lat);
    mooring->lat=lat;
    mooring->lon=lon;
    mooring->depth=0;
    mooring->name=strdup(name);
    }
  
/*----------------------------------------------------------------------------
  read data */
  int nrecords = mgr_line_count(filename)-1;

  double z;
  double *t        = new double[nrecords];
  double *sealevel = new double[nrecords];

  int k = 0;
  int year = 0, month = 0, day = 0, hour = 0, minute = 0;
  double second = 0;

  for(int i = 0; i < nrecords; i++){
    std::getline( input, line );
    if(line[0]=='#') continue;
    
    int cnt = scan_dmy_hms_z(line, &year, &month, &day, &hour, &minute, &second, &z, i);
    if( cnt != 7) continue;
    
    sealevel[k]=z/100.;
    t[k] = julian_day(month, day, year)  - CNES0jd + hour /24. + minute / (24. * 60.) + second / (24.*3600.);
//     if(i>0) {
//       double zz=(sealevel[k]-sealevel[k-1])/(t[k]-t[k-1])/24.;
//       if ( (fabs(zz) > 1000.0) && (sealevel[k-1]!=serie->mask) ) {
//         sealevel[k]=serie->mask;
//         }
//       }
    k++;
    }
//   for (int i=0; i< k; i++) {
//     int kmin=max(i-10,0);
//     int kmax=min(i+10,k);
//     statistic_t s=get_statistics(&(sealevel[kmin]),serie->mask,kmax-kmin-1,0);
//     if(sealevel[i] != serie->mask) {
//       if(abs(sealevel[i]-s.mean) > 3*s.std) {
//         sealevel[i]=serie->mask;
//         }
//       }
//     }
//
//   for (int i=0; i< k; i++) {
//     int kmin=max(i-10,0);
//     int kmax=min(i+10,k);
//     statistic_t s=get_statistics(&(sealevel[kmin]),serie->mask,kmax-kmin-1,0);
//     if(sealevel[i] != serie->mask) {
//       if(abs(sealevel[i]-s.mean) > 3*s.std) {
//         sealevel[i]=serie->mask;
//         }
//       }
//     }

  serie->x=new double*[1];
  serie->x[0]=sealevel;
//   serie->x[1]=0;
//   serie->x[2]=0;
  
  serie->t=t;
  serie->n=k;
  serie->nparam=1;
  serie->first=poctime_getdatecnes(t[0],'d');
  
// #if 0
//   for (int i=0; i< serie->n; i++) {
//     int kmin=max(i-50,0);
//     int kmax=min(i+50,serie->n);
// //    statistic_t s=get_statistics(&(residuals[kmin]),serie->mask,kmax-kmin-1,0);
//     statistic_t s=get_statistics(&(serie->x[2][kmin]),serie->mask,kmax-kmin-1,0);
//     if(serie->x[0][i] != serie->mask) {
// //      serie->x[2][i]=serie->x[0][i]-residuals[i]-mean;
//       if(abs(residuals[i]-s.mean) > 6*s.std) {
//         serie->x[0][i]=serie->mask;
//         }
//       }
//     }
// #endif
//   for (int i=0; i< serie->n; i++) {
//     if(serie->x[0][i] == serie->mask) {
//       serie->x[1][i]=serie->mask;
//       serie->x[2][i]=serie->mask;
//       }
//     }
    
  input.close();


  return(1);
}
