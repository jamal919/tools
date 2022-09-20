
/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/
#include <config.h>

#include <stdio.h>
#include <string.h>

#include "tools-structures.h"
#include "functions.h"
#include "poc-time.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int julian_day(const date_t &date)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=julian_day/*avoid perlReplace.pl*/(date.month,date.day,date.year);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int julian_day(int month, int day, int year)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  routine returns the julian day number which begins at noon
  of the calendar date specified by month mm, day id, and year iyyy.
  positive year signifies a.d.; negative, b.c.  (remember that the
  year after 1 b.c. was 1 a.d.
  
  this routine has been lifted directly from the book
  press et al., numerical recipes, cambridge univ. press, 1986.
  
  routine handles changeover to gregorian calendar on oct. 15, 1582.
  
  note: to get the corresponding modified julian day,
        set mjd = julday - 2400001.
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double jd;
  int gregorian=15+31*(10+12*1582);
  int jy,jm,ja, tmp_iyyy,out;

  tmp_iyyy=year;

  if (tmp_iyyy == 0)
    TRAP_ERR_EXIT(99,"Year 0 does not exist!\n");
/*   if (tmp_iyyy >= 2000) printf("Probleme avec le passage a l'annee 2000 \n"); */

  if (tmp_iyyy < 0) tmp_iyyy = tmp_iyyy + 1;
  if (month > 2) {
    jy = tmp_iyyy;
    jm = month + 1;
    }
  else {
    jy = tmp_iyyy - 1;
    jm = month + 13;
    }

  jd = floor(365.25*jy) + floor(30.6001*jm) + day + 1720995;

  if (day+31*(month+12*tmp_iyyy) >= gregorian) {
    ja = (int)( floor (0.01*jy) );
    jd = jd + 2 - ja +  floor(0.25*ja);
    }
  out=(int) jd;
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mjd(int month, int day, int year)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  routine returns the modified julian day number which begins at noon
  of the calendar date specified by month mm, day id, and year iyyy.
  positive year signifies a.d.; negative, b.c.  (remember that the
  year after 1 b.c. was 1 a.d.
  
  this routine has been lifted directly from the book
  press et al., numerical recipes, cambridge univ. press, 1986.
  
  routine handles changeover to gregorian calendar on oct. 15, 1582.
  
  note: to get the corresponding modified julian day,
        set mjd = julday - 2400001.
        
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int jd;

  jd=julian_day(month, day, year)- 2400001;
  return(jd);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int day_in_year(const date_t &date)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------

  compute time elapsed in year in integer days

!-----------------------------------------------------------------------*/

{
  int jd;
  date_t tmp=date_t(date.year,1,1,0);

  jd=julian_day(date)- julian_day(tmp)+1;

  return(jd);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ndays_in_year(int year)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *-----------------------------------------------------------------------

  compute time elapsed in year in integer days

!-----------------------------------------------------------------------*/

{
  int jd;

  jd=julian_day(1,1,year+1)- julian_day(1, 1,year);

  return(jd);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void calendary(int njd,date_t *actual)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nj;
  int na,nm,nd;

/* *-----------------------------------------------------------------------
! Input = njd: nombre ecoules depuis le 1er  Janvier 1950, 0 heure
! Output= nd,nm,na: jour, mois annee calendaire
!-----------------------------------------------------------------------*/

  nj=njd;

  na=0;
  if(nj<0) {
    while(nj<0) {
      na--;
      nj+=ndays_in_year(na+1950);
      }
    }
  else {
    while(nj>=0) {
      nj-=ndays_in_year(na+1950);
      na++;
      }
    na--;
    nj+=ndays_in_year(na+1950);
    }

  nj=nj+1;

  nm=1;
  while(nj>poctime_dpm(nm, na+1950)) {
    nj-=poctime_dpm(nm, na+1950);
    nm++;
    }
 
  nd=nj;

  actual->year=na+1950;
  actual->month=nm;
  actual->day=nd;
  actual->second=0;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void calendary_obsolete(int njd,date_t *actual)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int njul,nj,nb,nm1,nj3,m,j,ndj;
  int na,naj,nm,nd;
  int n[12]= {31,28,31,30,31,30,31,31,30,31,30,31};

/* *-----------------------------------------------------------------------
! Input = njd: nombre ecoules depuis le 1er  Janvier 1950, 0 heure
! Output= nd,nm,na: jour, mois annee calendaire
!-----------------------------------------------------------------------*/

  njul=njd+1;
  naj=njul/365;
/* *----------------------------------------------------------------------------
  patch for negative year */
  if(naj<0) naj--;
  nj=njul-naj*365;

/* *----------------------------------------------------------------------------
  patch for negative year */
  if(nj<0) {
    naj--;
    nj=njul-naj*365;
    }

  na=njul/365;
  if(na<0) na--;

  nb=(na+1)/ 4;
/* *----------------------------------------------------------------------------
  patch for negative year */
  if(na<0) {
    nb=(na-2)/ 4;
    }

  if(na+1950<1900) {
/* *----------------------------------------------------------------------------
    centuries are not leap year, except millenaries */
    nb++;
    }

/* *----------------------------------------------------------------------------
  compensate for leap years */
  nj=nj-nb;

/* *----------------------------------------------------------------------------
  patch for negative year */
  if(naj<0) {
    if(nj>ndays_in_year(naj+1950)) {
      nj-=ndays_in_year(naj+1950);
      naj++;
      }
    }

/* *----------------------------------------------------------------------------
  centuries extra patch */
  if(naj+1950==1900) {
    if((nj>13) && (nj<61)) {
      nj--;
      }
    }

  if (nj >  0) goto a1;

/* *----------------------------------------------------------------------------
  initial coding: simply dismiss case if negative time */
  naj=naj+1949;
  nm=12;
  nd=nj+31;
  goto a9000;

a1:
/* *----------------------------------------------------------------------------
  j is zero if bissextile year */
  j=naj-2-nb*4;
  naj=naj+1950;

  if (j < 0) goto a5000;

// a4000:
/* *----------------------------------------------------------------------------
  leap year, later than end of february, special initialization required */
  if (60-nj < 0)  goto  a4500;
/* *----------------------------------------------------------------------------
  leap year, last day fo february, date found*/
  if (60-nj == 0) goto  a7000;
/* *----------------------------------------------------------------------------
  leap year, earlier than end of february (so equivalent to common year) */
  goto a5000;

a4500:
/* *----------------------------------------------------------------------------
  leap year, february passed, initialize and continue treatment as a common year */
  nm1=60;

  m=3;
  goto a6000;

a5000:
/* *----------------------------------------------------------------------------
  common year, initialize loop for month computation */
  nm1=0;

  m=1;

a6000:
  ndj=nm1+n[m-1];

  nj3=nj-ndj;

  if (nj3 <= 0) goto a8000;

// a6500:
/* *----------------------------------------------------------------------------
  increment month */
  m=m+1;

  nm1=ndj;
  goto a6000;

a7000:
/* *----------------------------------------------------------------------------
  leap year, last day fo february*/
  nm=2;

  nd=29;

  goto a9000;

a8000:
  nm=m;

  nd=nj-nm1;

a9000:
  actual->year=naj;
  actual->month=nm;
  actual->day=nd;
  actual->second=0;

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poctime_dpm(int month, int year)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  dpm= days per months
!-----------------------------------------------------------------------*/
{
  int jd;

  if (month < 12)
    jd=julian_day(month+1, 1, year)- julian_day(month, 1, year);
  else
    jd=julian_day(1, 1, year+1)- julian_day(month, 1, year);
  return(jd);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// double cnes_stime(date_t actual)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int N;
//   float second;
//   double time;
//
// /*----------------------------------------------------------------------
//   return cnes time in second corresponding to calendar date*/
//
//   N=julian_day(actual) -julian_day(1,1,1950);
//
//   time=N*24.*3600.+actual.second;
//
//   return (time);
// }
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// double cnes_dtime(date_t actual)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int N;
//   float second;
//   double time;
//
// /*----------------------------------------------------------------------
//   return cnes time in day corresponding to calendar date*/
//
//   N=julian_day(actual) -julian_day(1,1,1950);
//
//   time=N+(actual.second)/(24*3600.);
//
//   return (time);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double time_elapsed(date_t date1, date_t date2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  /**
  \brief return in t the time elapsed between date1 and date2

  @param [in]     date1
  @param [in]     date2
  @return double     the time elapsed between date1 and date2
  */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

{
  int N;
  double t;

  N = julian_day(date2.month, date2.day, date2.year)
        - julian_day(date1.month, date1.day, date1.year);

  t = N * 24. * 3600. + date2.second - date1.second;
  return (t);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double cnes_time(date_t actual, const char units)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

get CNES time in requested unit of given date

  actual : given date. An invalid date (see isnad()) will give an error.
  unit : requested unit : 's' for seconds,'h' for hours,'d' for days or 'y' for Julian years. Anything else will give an error.

returns CNES time in requested unit of given date or NAN if error

For the revert, see poctime_getdatecnes().

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int N;
  double time;
  double factor2;

  switch (units) {
    case 's':
    case 'S':
      factor2 = 1.;
      break;
    case 'h':
    case 'H':
      factor2 = 3600.;
      break;
    case 'd':
    case 'D':
      factor2 = 3600. *24.;
      break;
    case 'y':
    case 'Y':
      factor2 = 3600. *24. *365.25;
      break;
    default:
      TRAP_ERR_EXIT(-1, "invalid time unit\n");
      return NAN;
    }

  if(actual.year==0)
    return NAN;
  ///\note CNES time origin is 1950-01-01 00:00:00
  N=julian_day(actual)-julian_day(1,1,1950);

  time=N * 24. * 3600. +actual.second;

  time/=factor2;

  return (time);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poctime_convert(double *t, int n, const char units_in, const char units_out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k;
  double factor, factor1, factor2;
  
  if(units_in==units_out) return 0;

/*------------------------------------------------------------------------------
  return cnes time in requested unit corresponding to calendar date*/
  switch (units_in) {
    case 's':
    case 'S':
      factor1 = 1.;
      break;
    case 'h':
    case 'H':
      factor1 = 3600.;
      break;
    case 'd':
    case 'D':
      factor1 = 24. *3600.;
      break;
    case 'y':
    case 'Y':
      factor1 = 24. *3600. *365.25;
      break;
    default:
      return(-1);
      break;
    }
  switch (units_out) {
    case 's':
    case 'S':
      factor2 = 1.;
      break;
    case 'h':
    case 'H':
      factor2 = 3600.;
      break;
    case 'd':
    case 'D':
      factor2 = 24. *3600.;
      break;
    case 'y':
    case 'Y':
      factor2 = 24. *3600. *365.25;
      break;
    default:
      return(-1);
      break;
    }
  factor=factor2/factor1;
  for (k=0;k<n;k++) t[k]/=factor;
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

double ellapsed_time(date_t start, date_t final, const char units)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double time;
  time=cnes_time(final, units) - cnes_time(start, units);
  return (time);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// void getcnesdate(double t,date_t *actual)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int N,nday;
//   float second;
//
// /* *----------------------------------------------------------------------
//   Get the calendary date corresponding to cnes time t (in decimal hours
//   elapsed from the 1/1/1950). The residual in seconds is rounded
// ----------------------------------------------------------------------*/
//
// //  nday=(int)( t/24. );
//   nday=(int)( floor(t/24.) );
//
//   calendary(nday,actual);
//   second=(t-nday*24.)*3600.;
//
//   actual->second=floor(second+.5);
//
// }

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// date_t cnesdate(double t, date_t *date, char fmt)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
// /* *----------------------------------------------------------------------
//   Get the calendary date corresponding to cnes time t (
//   elapsed from the 1/1/1950). The residual in seconds is rounded
// ----------------------------------------------------------------------*/
// {
//   int N,nday;
//   float second;
//   date_t actual;
//
//   double factor1,factor2;
//
//   switch (fmt) {
//     case 's':
//     case 'S':
//       factor1=24.*3600.;
//       factor2=1.;
//       break;
//     case 'h':
//     case 'H':
//       factor1=24.;
//       factor2=3600.;
//       break;
//     case 'd':
//     case 'D':
//       factor1=1.;
//       factor2=24.*3600.;
//       break;
//     default:
//       TRAP_ERR_EXIT(-1,"exiting\n");
//     }
//
//   if(t>0) {
//     nday=(int)( floor(t/factor1) );
//     second=(t-nday*factor1)*factor2;
//     }
//   else {
//     nday=(int)( floor(t/factor1) );
//     second=(t-nday*factor1)*factor2;
//     }
//
//   calendary(nday,&actual);
//
//   actual.second=floor(second+.5);
//
//   *date=actual;
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  date_t poctime_getdatecnes(double t, char unit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

convert double since 1/1/1950 to date_t

  t : s since 1 Jan 1950
  unit : 's', 'h' or 'd'
  
For the revert, see cnes_time().

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int nday;
  float second;
  date_t actual;

  double factor1=NAN,factor2=NAN;
  
  if(isnan(t))
    return NADate;

  switch (unit) {
    case 's':
    case 'S':
      factor1=24.*3600.;
      factor2=1.;
      break;
    case 'h':
    case 'H':
      factor1=24.;
      factor2=3600.;
      break;
    case 'd':
    case 'D':
      factor1=1.;
      factor2=24.*3600.;
      break;
    default:
      TRAP_ERR_EXIT(-1,"exiting\n");
    }

  nday=(int)( floor(t/factor1) );
  second=(t-nday*factor1)*factor2;

  calendary(nday,&actual);

  actual.second=second;

  return(actual);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  date_t poctime_getdatenasa(double t, char unit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nday;
  float second;
  date_t actual;
/*------------------------------------------------------------------------------
  Get the calendary date corresponding to nasa time t (in decimal hours
  elapsed from the 1/1/1985). The residual in seconds is rounded
----------------------------------------------------------------------*/
  switch (unit)
    {
    case 's':
    case 'S':
      nday=(int)( t/24./3600. );
      break;
    case 'h':
    case 'H':
      nday=(int)( t/24. );
      break;
    case 'd':
    case 'D':
      nday=(int)(t);
      break;
    default:
      TRAP_ERR_EXIT(-1,"exiting\n");
    }

  nday=(int)( floor(t/24.)+mjd(1,1,1985)-mjd(1,1,1950) );

  calendary(nday,&actual);
  second=(t-floor(t/24.)*24.)*3600.;

  actual.second=second;

  return(actual);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  date_t poctime_getdate(double t, date_t reference, char unit, double *cnestime)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  Get the calendary date corresponding to  time t (in format unit
  elapsed from reference). The residual in seconds is rounded
----------------------------------------------------------------------*/
{
  int nday;
  float second;
  double factor1=NAN,factor2=NAN;
  date_t actual;

  switch (unit)
    {
    case 's':
    case 'S':
      factor1= 24.* 3600.;
      factor2= 1.;
      break;
    case 'h':
    case 'H':
      factor1= 24.;
      factor2=3600.;
      break;
    case 'd':
    case 'D':
      factor1= 1.;
      factor2= 24.* 3600.;
      break;
    default:
      TRAP_ERR_EXIT(-1,"exiting\n");
    }

  t+=reference.second;
  nday=(int)( floor(t/factor1) );
  second=(t-nday*factor1)*factor2;

  nday=nday+mjd(reference.month,reference.day,reference.year)-mjd(1,1,1950);

  calendary(nday,&actual);

  actual.second=second;

  *cnestime=nday+second/( 24.* 3600.);

  return(actual);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *sgetcnesdate(double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  date_t actual;
  char *s=NULL;

/*------------------------------------------------------------------------------
  t is cnes time in hours
----------------------------------------------------------------------*/

  actual=poctime_getdatecnes(t,'h');
  
  s=sgetdate(actual);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *poctime_sdate_cnes(double t, char unit, const char sep)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///CNES date to string as "yyyy/mm/dd HH:MM:SS"
/**
\param t CNES date
\param unit
\param separator character to separate date and time. Default: space
\returns the string allocated with \c new
*/
/*----------------------------------------------------------------------------*/
{
  date_t actual;
  char *s=NULL;
  
  if(t==0. && unit=='\0')
    unit='s';//whatever
  
  exitIfNull(s=new char[64]);
  if(isnan(t) || (unit=='\0' && t!=0.)){
    sprintf(s,"NADate");
    return s;
    }
  if(not isfinite(t)){
    sprintf(s,"%g",t);
    return s;
    }
  
  actual=poctime_getdatecnes(t,unit);

  s=sgetdate(actual,sep);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *sgetdate(date_t actual, const char sep)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int hour=-1,minute;
  double second;
  char *s=NULL;
  
  fct_hms(actual.second/36e2,&hour,&minute,&second);

  exitIfNull(
    s=new char[64]
    );

  //note: float seconds of day is precise to 0.0078
  sprintf(s,"%04d/%02d/%02d%c%02d:%02d:%04.1f",actual.year,actual.month,actual.day,sep,hour,minute,second);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *sgetnasadate(double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  date_t actual;
  char *s=NULL;

/*------------------------------------------------------------------------------
  Get the calendary date corresponding to absolute time t_reference+t
----------------------------------------------------------------------*/
  actual=poctime_getdatenasa(t,'s');
  
  s=sgetdate(actual);

  return(s);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

char *sgetnasadatef_01(float *t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double tt;
  tt=*t;
  return(sgetnasadate(tt));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void sgetnasadatef_02(float *t, char *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double tt;
  tt=*t*24.0;
  strcpy(s,sgetnasadate(tt));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool poctime_date_valid(date_t *date)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if (date->year < 1500 || date->year > 2500 ||
      date->day > 31 || date->day < 1 ||
      date->month > 12 || date->month < 1)
    return false;
  return true;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  date_t poctime_scan_date(const char *date_str, date_t *date_default)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* ----------------------------------------------------------------------

Convert date strings to date_t. Following formats are accepted:
 yyyy/mm/dd
 dd/mm/yyyy

'/', dd and mm are optional
   ---------------------------------------------------------------------- */
{
  date_t tmp;
  int nitems;

  /* just a year */
  if (strlen(date_str) == 4 && isdigit(date_str[0]) && isdigit(date_str[1]) && isdigit(date_str[2]) && isdigit(date_str[3]) ) {
    tmp=*date_default;
    nitems=sscanf(date_str,"%d", &tmp.year);
    if(nitems==1 && poctime_date_valid(&tmp))
      return tmp;
  }
  
  /* there is a '/' */
  if (strchr(date_str,'/')) {
    /* read format dd/mm/yyyy */
    tmp=*date_default;
    nitems=sscanf(date_str,"%2d/%2d/%4d",&tmp.day,&tmp.month,&tmp.year);
    if (nitems > 1 && poctime_date_valid(&tmp))
      return tmp;
    
    /* try yyyy/mm/dd */
    tmp=*date_default;
    nitems=sscanf(date_str,"%4d/%2d/%2d",&tmp.year,&tmp.month,&tmp.day);
    if (nitems > 1 && poctime_date_valid(&tmp))
      return tmp;
  }

  /* no '/' */
  else {
    /* read format ddmmyyyy */
    tmp=*date_default;
    nitems=sscanf(date_str,"%2d%2d%4d",&tmp.day,&tmp.month,&tmp.year);
    if (nitems > 1 && poctime_date_valid(&tmp))
      return tmp;
    
    /* try yyyymmdd */
    tmp=*date_default;
    nitems=sscanf(date_str,"%4d%2d%2d",&tmp.year,&tmp.month,&tmp.day);
    if (nitems > 1 && poctime_date_valid(&tmp))
      return tmp;
  }

  STDERR_BASE_LINE("wrong date format %s", date_str);
  return NADate;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_poctime_scan_date_help(int end_date_mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// print user help for poctime_scan_date()
/*----------------------------------------------------------------------------*/
{
  printf("\n"
    "DATE FORMATS\n"
    "The following formats are accepted:\n"
    "  yyyy/mm/dd HH:MM:SS.SSS\n"
    "  dd/mm/yyyy HH:MM:SS.SSS\n"
    "  mm/yyyy\n"
    "The least significant part (second, minute, hour, day, month) is always optional.\n"
    "The separator character can be any non-number character (anything other than -+0123456789.eEfF)\n"
    );
  
  if(end_date_mode>=0){
    printf(
      "End dates are inclusive. This means that end dates specified as:\n"
      "  `2000`             will be interpreted as `2001/01/01 00:00:00`\n"
      "  `2000/01`          will be interpreted as `2000/02/01 00:00:00`\n"
      "  `2000/01/01`       will be interpreted as `2000/01/02 00:00:00`\n"
      "  `2000/01/01 00`    will be interpreted as `2000/01/01 01:00:00`\n"
      "  `2000/01/01 00:00` will be interpreted as `2000/01/01 00:01:00`\n"
      );
    }
  
  if(end_date_mode>=1){
    printf(
      "End dates default to the start date in an inclusive way.\n"
      );
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int poctime_scan_date(const string &str, date_t *start, date_t *end, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// convert date string to start and end date_t
/**
\param str
\param *start
\param *end inclusive. If already set, leave unchanged
\return 0 if the format is OK, -1 otherwise.

\sa print_poctime_scan_date_help().

Replace with
\code
sed -i~ -re 's/sscanf\(([^,]+?),"%d\/%d\/%d",&([^,]+?)\.day,&\2\.month,&\2\.year\);/poctime_scan_date(\1,\&\2,0);/' *.cpp
\endcode
*/
/*----------------------------------------------------------------------------*/
{
  int i,nitems,count[6];
  int hour,minute,jd;
  date_t *date;
  float second;
  const char *s=str.c_str();
  
  for(i=0;i<2;i++){
    
    /* select output */
    if(i==0){
      date=start;
      }
    else if(i==1){
      if(start!=0 and end!=0 and isad(*end)) /* If already set */
        /* leave unchanged */
        continue;
      date=end;
      }
    else TRAP_ERR_EXIT(ENOEXEC,"coding error\n");
    
    if(date==0)
      continue;
    
    /* initialise */
    *date=date_t(0,1,1,0.0);
    hour=minute=0;
    second=0.f;
    
    /* scan input */
    aset(count,6,0);
    sscanf(s,"%d%n%*c%d%n%*c%d%n%*c%d%n%*c%d%n%*c%g%n",&date->year,&count[0],&date->month,&count[1],&date->day,&count[2],&hour,&count[3],&minute,&count[4],&second,&count[5]);
    
    if(count[0]<=0){
      *date=NADate;
      TRAP_ERR_RETURN(-1,verbose,"could not scan date %s\n",s);
      }
    
    /* remove count offsets */
    for(nitems=5;nitems>0;nitems--)
      if(nitems>0)
        count[nitems]-=count[nitems-1]+1;
    
    /* count scanned fields */
    for(nitems=0;nitems<6 and count[nitems]>0;nitems++);
    
    if(nitems==2 and count[0]<count[1])/* mm/yyyy format */
      swapValues(&date->year,&date->month);
    
    if(nitems>=3 and count[0]<count[2])/* dd/mm/yyyy format */
      swapValues(&date->year,&date->day);
    
    /* final processings */
    second+=(hour*60.+minute)*60.;
    jd=julian_day(*date);
    
    if(i>0){
      switch(nitems){
      case 1:
        /* year */
        jd+=ndays_in_year(date->year);
        break;
      case 2:
        /* month */
        jd+=poctime_dpm(date->month,date->year);
        break;
      case 3:
        /* day */
        jd+=1;
        break;
      case 4:
        /* hour */
        second+=3600.f;
        break;
      case 5:
        /* minute */
        second+=60.f;
        break;
        }
      }
    
    calendary(jd-CNES0jd,date);
    date->second=second;
    
    if(verbose>0)
      STDERR_BASE_LINE("[%d]%s\n",i,sgetdate(*date));
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

date_t poctime_scan_start_date(const char *date_str)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* ----------------------------------------------------------------------

Convert date strings to date_t. Following formats are accepted:
 yyyy/mm/dd
 dd/mm/yyyy

dd and mm are optional, and replaced by 01
   ---------------------------------------------------------------------- */
{
  date_t date=date_t(0,1,1,0.0);
  return poctime_scan_date(date_str, &date);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  date_t poctime_scan_final_date(const char *date_str)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* ----------------------------------------------------------------------

Convert date strings to date_t. Following formats are accepted:
 yyyy/mm/dd
 dd/mm/yyyy

mm is optional and replaced by 12
dd is optional and replaced by the last day of the month

   ---------------------------------------------------------------------- */
{
  date_t date=date_t(0,12,31,0.);
  return poctime_scan_date(date_str, &date);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  date_t poctime_add_step(date_t date, date_t step, date_t stop, date_t *date_end)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  return staring date of next step
     date_end : date at end of the step (ie, 1 day before date of next step)
     
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  date_t date_next(date);
  int d_date;

  *date_end = date;

/*------------------------------------------------------------------------------
  add years */
  if (step.year!=0) {

    // date for next step
    date_next.year += step.year;

    // date of end of step : convert to number of days
    d_date = cnes_time(date_next, 'd');
    calendary(d_date, date_end);
    }

/*------------------------------------------------------------------------------
  add months */
  if (step.month!=0) {

    // date for next step
    date_next.month += step.month;
    while (date_next.month >12) {
      date_next.month -= 12;
      date_next.year += 1;
      }

    // date of end of step : convert to number of days
    d_date = cnes_time(date_next, 'd');
    calendary(d_date, date_end);
    }
  
/*------------------------------------------------------------------------------
  there is a limit */
  if (isad(stop)) {
    *date_end = poctime_soonest(*date_end, stop);
    if (poctime_is_after(date_next, stop))
      return NADate;
    }

  return date_next;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int poctime_count_steps(date_t date_start, date_t date_final, date_t step)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  date_t next = date_start;
  date_t end;
  int nb_steps=1;

  if (date_start > date_final)
    return 0;
  while( isad(next = poctime_add_step(next, step, date_final, &end)) )
    nb_steps++;

#if DEBUG
  printf ("poctime_count_steps nb_steps=%d\n", nb_steps);
#endif
  return nb_steps;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

date_t poctime_soonest(date_t date0, date_t date1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int d0, d1;
//   d0 = julian_day(date0.month, date0.day, date0.year);
//   d1 = julian_day(date1.month, date1.day, date1.year);
//   if (d0 <= d1)
//     return date0;
//   return date1;
  if (date1 > date0)
    return date0;
  return date1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

date_t poctime_latest(date_t date0, date_t date1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int d0, d1;
//   d0 = julian_day(date0.month, date0.day, date0.year);
//   d1 = julian_day(date1.month, date1.day, date1.year);
//   if (d0 < d1)
//     return date1;
//   return date0;
  if (date0 > date1)
    return date0;
  return date1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

bool poctime_is_equal(date_t date0, date_t date1)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   if (date0.year == date1.year &&
//       date0.month == date1.month &&
//       date0.day == date1.day)
//     return true;
//   return false;
  return date0 == date1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

bool poctime_is_before(date_t before, date_t after)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int d0, d1;
//   d0 = julian_day(before.month, before.day, before.year);
//   d1 = julian_day(after.month, after.day, after.year);
//   if (d0 < d1)
//     return true;
  return before < after;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

bool poctime_is_after(date_t after, date_t before)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
//   int d0, d1;
//   d0 = julian_day(before.month, before.day, before.year);
//   d1 = julian_day(after.month, after.day, after.year);
//   if (d0 > d1)
//     return true;
//   return false;
  return after > before;

}
