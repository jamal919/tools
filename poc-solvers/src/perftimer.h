
/**************************************************************************

  POC-SOLVERS interface, 2006-2019

  Part of the SIROCCO national service (INSU, France)

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France

***************************************************************************/

#include <string>

#include "mpi-reducers.h"
#include "solvers-statistics.h"

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  acquire_start and acquire_stop must be implemented just before and just after
  the tested function call

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

class perftimer_t {
  string baseLineFunc;
  struct timeval t;
  double *durations;
  double start, stop;
  
public:
  
  string report;
  statistic_t<double> global;
  int id,ntics,tic;
  
public:
  
  void init() {
    durations=0;
    ntics=0; tic=0;
    start=NAN;
    stop=NAN;
    }
    
  void init(int range) {
    ntics=range;
    durations=new double[ntics];
    }
    
  perftimer_t() {
    init();
    }
    
  void process() {
    statistic_t<double> s;
    s=compute_statistics(durations, NAN, ntics);
    global.mean=P_average(s.mean); 
    global.max=P_max(s.max);
    char *message=new char[1024];
    sprintf(message, "average=%lfms max=%lfms (ntics=%d)", 1000*global.mean, 1000*global.max, ntics);
    report=message;
    delete[] message;
    }
    
//   void print() {
//     statistic_t<double> s, global;
//     s.compute_statistics(durations, NAN, ntics);
//     global.mean=P_average(s.mean); 
//     global.max=P_max(s.max);
// //     printf("%s : solver %d (neq=%d) solving time average=%lfms max=%lfms \n", __func__, solver->id, 0, average, max);
//     }
    
  void acquire_start() {
    double t;
    if(ntics<1) return;
    stop=NAN;
#ifdef HAVE_MPI
    start = MPI_Wtime();
#endif
    if(tic==ntics) tic=0;
    }
    
  bool acquire_stop() {
    double t;
    bool ready=false;
    if(ntics<1) return(ready);
#ifdef HAVE_MPI
    stop = MPI_Wtime();
#endif
    t=stop-start;
    durations[tic]=t;
    start=NAN;
    stop=NAN;
    tic++;
    if(tic==ntics) {
      process();
      ready=true;
      }
    return(ready);
    }
    
  void check() {
    if(tic>=ntics) {
      tic=0;
      }
    }
    
  void destroy() {
    ptr_delete(durations);
    init();
    }
    
  ~perftimer_t(){
    destroy();
    };
  };


