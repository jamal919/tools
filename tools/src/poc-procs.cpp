
/*******************************************************************************

  T-UGO tools, 2006-2016

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief definition of usefull functions missing from standard libraries: process and optimisation
*/
/*----------------------------------------------------------------------------*/

#include <stdio.h> /* for printf() */
#include <stdlib.h> /* for system() */
#include <sys/time.h> /* for gettimeofday() */
#include <time.h> /* for gmtime_r */
#include <unistd.h> /* for usleep() */
#include <limits.h> /* for INT_MAX */
#include <fstream>

#include "../config.h" /* for ERROR_ON_EXPIRY */
#include "functions.h" /* safely includes omp.h */
#include "poc-time-classes.h" /* for date_t */
#include "poc-assertions.h"
#include "poc-array.hpp" /* for deletep() */


time_t lastTime=time(NULL);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  time_t timeIsOld()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  time_t thisTime=time(NULL);
  
  if(thisTime>lastTime){
    #pragma omp critical(timeIsOldIsThreadSafe)
    {
    if(thisTime>lastTime)
      lastTime=thisTime;
    else
      thisTime=0;
    }
    
    return thisTime;
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double difftime(const struct timeval &b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
///difftime overload: difference between current time and argument
/**
\code
struct timeval before;
gettimeofday(&before);
//... some I/O code ...
iot+=difftime(before);
\endcode
*/
{
  struct timeval a;
  gettimeofday(&a);
  return (double)(a.tv_sec-b.tv_sec)+1e-6*(a.tv_usec-b.tv_usec);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int isnad(const date_t &date)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// is Not A Date
/**
\note As negative months and days are tolerated by POSIX date functions,
it tests the values of neither date_t::day nor date_t::month .
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  return isnan(date.second)||(date.year==0); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int isad(const date_t &date)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// is A Date
/**
This is faily permissive : it returns the oposite of isnad(). See its note.
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  return !isnad(date); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  date_t timeval2date(const struct timeval & tv)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  struct tm b;
  gmtime_r(&tv.tv_sec,&b);
  
  const float S=(b.tm_hour*60.f+b.tm_min)*60.f+b.tm_sec+tv.tv_usec*1e-6f;
  
  const date_t date(b.tm_year+1900,b.tm_mon+1,b.tm_mday,S);
  
  return date;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ymd2ymd(int ymd,int *y,int *m,int *d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const int ym=ymd/100;
  *d=ymd%100;
  *m=ym%100;
  *y=ym/100;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool expire_(const char *file,int line,int ymde, int ymdr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// make sure out-dated code is not executed unknowingly
/**
Do no call directly. Call #expire instead.

\param *file set to \c __FILE__
\param line  set to \c __LINE__
\param ymde date of expiry  as YYYYMMDD
\param ymdr date of removal as YYYYMMDD. Will exit if past this date.
\returns whether past date of expiry
*/
/*----------------------------------------------------------------------------*/
{
  struct timeval now;
  gettimeofday(&now);
  const date_t now_date=timeval2date(now);
  
  int y,m,d;
  
  if(ymdr>0){
    ymd2ymd(ymdr,&y,&m,&d);
    const date_t removal(y,m,d);
    
    if(removal<now_date){
      fprintf(stderr,"*** %s:%d:SHOULD BE REMOVED SINCE %04d-%02d-%02d ***\n",file,line,y,m,d);
#if ERROR_ON_EXPIRY
      wexit(ENOEXEC);
#else
      return true;
#endif
      }
    }
  
  ymd2ymd(ymde,&y,&m,&d);
  const date_t expiry(y,m,d);
  
  if(now_date<expiry){
    fprintf(stderr,"*** %s:%d:WILL EXPIRE ON %04d-%02d-%02d UNLESS YOU DO SOMETHING ABOUT IT ***\n",file,line,y,m,d);
    return false;
    }
  
  fprintf(stderr,"\n*** PREVIOUS BEHAVIOUR %s:%d:EXPIRED ON %04d-%02d-%02d ***\n\n",file,line,y,m,d);
  return true;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cmp(const struct timeval &a,const struct timeval &b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  
  r=cmp(a.tv_sec,b.tv_sec);
  
  if(r==0)
    r=cmp(a.tv_usec,b.tv_usec);
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double difftime(struct timeval *b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///difftime overload: difference between current and argument
/**
\param[in,out] *b pointer to last time value, will be updated to current time
\code
gettimeofday(&before);
//... some I/O code ...
iot+=difftime(before);
\endcode
*/
/*----------------------------------------------------------------------------*/
{
  double difference;
  struct timeval a;
  gettimeofday(&a);
  difference=(double)(a.tv_sec-b->tv_sec)+1e-6*(a.tv_usec-b->tv_usec);
  *b=a;
  return difference;
}


#ifdef CLOCK_REALTIME
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double difftime(const struct timespec b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///difftime overload, that takes current time as argument
/**
\bug 2011-09-21 Damien ALLAIN : I am not so sure about the reliability of clock_gettime
\code
clock_gettime(CLOCK_REALTIME,&before);
//... some I/O code ...
iot+=difftime(before);
\endcode
*/
/*----------------------------------------------------------------------------*/
{
  struct timespec a;
  clock_gettime(CLOCK_REALTIME,&a);
  return (a.tv_sec-b.tv_sec)+1e-9*(a.tv_nsec-b.tv_nsec);
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

difftimer_t::difftimer_t(const char *file,int line,const char *func,int n0,int id0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ostringstream oss;
  oss<<strrchr0(file,'/')<<":"<<line<<" in "<<func;
  baseLineFunc=oss.str();
  
  gettimeofday(&t);
  
  id=id0;
  
  n=n0;
  j=0;
  
  if(n>0)
    durations=aset(n,0.);
  else
    durations=0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void difftimer_t::operator()()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(durations==0)return;
  
  if(j<n)
    durations[j]+=difftime(&t);
  j++;
}


int verbose_difftimer=1;
char difftimer_unit_prefix='\0';


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void difftimer_t::destroy()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(durations==0)return;
  
  operator()();
  
/*------------------------------------------------------------------------------
  help correcting coding error */
  if(j>n){
    printf(""+baseLineFunc+": replace %d by at least %d\n",n,j);
    wexit(ENOEXEC);
    }
  
/*------------------------------------------------------------------------------
  output */
  if(verbose_difftimer>0){
    int i;
    double sum=0;
    double factor=NAN;
    
    switch(difftimer_unit_prefix){
    case '\0':
      factor=1.;
      break;
    case 'm':
      factor=.001;
      break;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"difftimer_unit_prefix=%c not recognised\n",difftimer_unit_prefix);
    }
    
    for(i=0;i<j;i++){
      durations[i]/=factor;
      sum+=durations[i];
      }
    
    /* MAKE SURE YOU UPDATE diffTimer.py:getTimerStats
      WHEN YOU MODIFY THE LINE BELOW!!! */
    printf(""+baseLineFunc+":difftimer_t(,%d) %g ",id,sum);
    if(difftimer_unit_prefix) putchar(difftimer_unit_prefix);
    printf("s =");
    
    for(i=0;i<j;i++)
      printf(" %g",durations[i]);
    
    printf("\n");
    fflush(stdout);
    }
  
/*------------------------------------------------------------------------------
  clean-up */
  deletep(&durations);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int cmp(const struct timespec &a,const struct timespec &b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int r;
  
  r=cmp(a.tv_sec,b.tv_sec);
  
  if(r==0)
    r=cmp(a.tv_nsec,b.tv_nsec);
  
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialize_OPENMP(int requested, FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nprocs;
  int gOPENMP_nCPUs;
  int gOPENMP;

#ifdef OMP_H
  const char *nprocsEnv=getenv("OMP_NUM_THREADS");
  if(nprocsEnv==0){
    /* get maximum number of cores on host */
    nprocs=omp_get_num_procs();
    }
  else{
    sscanf(nprocsEnv,"%d",&nprocs);
    }
  
  if(requested==-1)
    requested=nprocs;
  
  gOPENMP_nCPUs=min(nprocs,requested);
  
  if(gOPENMP_nCPUs>1){
    omp_set_schedule(omp_sched_dynamic,-1);
    omp_set_num_threads(gOPENMP_nCPUs);
    gOPENMP=1;
    }
  else {
    gOPENMP=0;
    gOPENMP_nCPUs=1;
    omp_set_num_threads(gOPENMP_nCPUs);
    }
#else
  gOPENMP=0;
  gOPENMP_nCPUs=1;
#endif
  
  if(out!=0) {
    fprintf(out,"OPENMP activated=%d, %d requested, %d available, %d used\n",gOPENMP,requested,nprocs,gOPENMP_nCPUs);
    }
  return(gOPENMP_nCPUs);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialize_OPENMP(int requested, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int gOPENMP_nCPUs;
  FILE *out;
  
  if(verbose)
    out=stdout;
  else
    out=0;
  
  gOPENMP_nCPUs=initialize_OPENMP(requested,out);
  
  return(gOPENMP_nCPUs);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int terminate_OPENMP(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  
#ifdef OMP_H
  omp_set_dynamic(0);
  omp_set_num_threads(1);
#endif
  
  printf("OPENMP terminated\n");
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_OPENMP_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints OpenMP help for a programme help function.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "ENVIRONMENT\n");
#if _OPENMP
  printf(
    "  This uses OpenMP version %d. Check the relevant API for more environment variables.\n",_OPENMP);
  printf(
    "  For information, OpenMP version 200505 and above are sensitive to the following variables :\n"
    "    OMP_SCHEDULE for the runtime schedule type and chunk size.\n"
    "    OMP_NUM_THREADS for the number of threads to use.\n"
    "    OMP_DYNAMIC for the dynamic adjustment of threads to use.\n"
    "    OMP_NESTED to enable or disable nested parallelism.\n"
    "  For example, in bash :\n"
    "    OMP_NUM_THREADS=6 %s ...\n",prog_name);
  printf("\n"
    "  If you are running on a machine with already loaded CPUs, you SHOULD take only a number of CPUs equivalent to the number that will remain free, otherwise the program may grind to an equivalent halt as soon as it parallelises.\n"
    );
#else
  printf(
    "  This version is not compiled with OpenMP support.\n");
#endif
  /** \endcode */
}


ostringstream assertionCmd;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int launchDebugCmd_(const char *file,int line)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string assertionCmdStr=assertionCmd.str();
  int status=-1;
  
  STDERR_BASE_LINE_FUNC("(from %s:%d) Running "+assertionCmdStr+" ======\n",file,line);
  if(assertionCmdStr.length()>0){
    status=system(assertionCmdStr.c_str());
    STDERR_BASE_LINE("returned %d ==================================\n",status);
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void get_cpus_stats(cpus_stats_t * stats)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ifstream f("/proc/stat");
  string s;
  int i,j;
  int64_t x;
  
  /* skip first line */
  f.ignore(INT_MAX,'\n');
  
  for(i=0;;i++){
    getline(f,s);
    
    if(strncmp("cpu",s)!=0)
      break;
    
    if(i>=stats->size())
      stats->resize(i+1);
    
    cpustats_t *stat=&(*stats)[i];
    stat->on=0;
    stat->idle=0;
    
    istringstream l(s);
    l>>s;
    for(j=0;;j++){
      l>>x;
      if(l.eof())break;
      
      /* search "/proc/stat" in man:/proc(5) */
      if(j==3)
        stat->idle+=x;
      else
        stat->on+=x;
      }
    }
  
  if(i<stats->size())
    stats->resize(i+1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fsleep(double t)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  useconds_t usec=t*1e6;
  usleep(usec);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int diff_cpus_stats(const cpus_stats_t & stats,vector<double> *ratios,FILE *out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  const int n=stats.size();
  cpus_stats_t newstats;
  get_cpus_stats(&newstats);
  
  if(n!=newstats.size())TRAP_ERR_RETURN(-1,1,"different number of processes : %d, then %d.\n",stats.size(),newstats.size());
  
  if(n!=ratios->size())
    ratios->resize(n);
  
  for(i=0;i<n;i++){
    const cpustats_t *stat0=&stats[i],*stat1=&newstats[i];
    const int64_t
      on=stat1->on-stat0->on,
      idle=stat1->idle-stat0->idle;
    
    double *ratio=&(*ratios)[i];
    
    *ratio=(double)on/(on+idle);
    
    if(out!=0){
      if(i>0)fprintf(out," ");
      fprintf(out,"cpu%d:%g%%",i,round(100.*(*ratio)));
      }
    
    }
  
  if(out!=0)
    fprintf(out,"\n");
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int diff_cpus_stats(const cpus_stats_t & stats,vector<double> *ratios,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  FILE *out;
  
  if(verbose)
    out=stdout;
  else
    out=0;
  
  result=diff_cpus_stats(stats,ratios,out);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int copyTo_tmpdir(string *path,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const char *tmpdir=getenv("TMPDIR");
  if(tmpdir==0) TRAP_ERR_RETURN(0,verbose,"%s: TMPDIR not set\n",__func__);
  
  
  string cmd="scp ";
  if(verbose>1)
    cmd+="-v ";
  cmd+=*path+" "+tmpdir;
  
  int status;
  if(verbose>0) STDERR_BASE_LINE("+ "+cmd+"\n");
  status=system(cmd.c_str());
  
  string newPath=tmpdir;
  if(strrncmp(newPath,"/")!=0)
    newPath+='/';
  newPath+=strrchr0(path->c_str());
  
  *path=newPath;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int copyTo_tmpdir(char **path,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  string s;
  
  s=*path;
  
  status=copyTo_tmpdir(&s,verbose);
  
  delete[]*path;
  
  *path=poc_strdup(s.c_str());
  
  return status;
}
