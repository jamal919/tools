
/**************************************************************************

  POC-SOLVERS interface, 2006-2012

  Part of the Unstructured Ocean Grid initiative and SYMPHONIE suite

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

***************************************************************************/

#define VERBOSE

#include <string>
#include <omp.h>
#include <complex>
#include <vector>

#include "solvers-interface.h"
#include "poc-solvers.h"


vector<SOLVER_REG_t <double> > dynamicsD;
vector<SOLVER_REG_t < complex<double> > > dynamicsZ;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int solver_register(int id, SOLVER_t<T> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
//   char *solver_name=get_SolverName(id);
//   solver.name = strdup(solver_name);
  
  switch(id) {
    case SOLVER_ID_DIAGONAL:     /* no real solver used for diagonal matrices */
      break;

    case SOLVER_ID_DOMESTIC:
      break;

    case SOLVER_ID_LAPACK:
      status=lapack_terminate(solver, verbose, debug);
      break;

    case SOLVER_ID_SUNPERF:
      break;

    case SOLVER_ID_UMFPACK:  
      status=umfpack_terminate(solver, verbose, debug);
      break;

    case SOLVER_ID_PETSC_GMRES:
      status=gmres_terminate(solver, verbose, debug);
      break;

    case SOLVER_ID_MUMPS: 
      status=mumps_terminate(solver, verbose, debug);
      break;

    case SOLVER_ID_MUMPS_SYM:
      status=mumps_terminate(solver, verbose, debug);
      break;

    case SOLVER_ID_SpDOMESTIC:
      break;

   case SOLVER_ID_HYPRE:
     break;
    
    case SOLVER_ID_HIPS:
      status=hips_terminate(solver, verbose, debug);
      break;  

    case SOLVER_ID_PASTIX:
      status=pastix_terminate(solver, verbose, debug);
      break;  

    default:
      return(-1);
      break;
    }


  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_identify_dynamic(SOLVER_t<double> & solver, int id, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  for(int k=0;k<dynamicsD.size();k++) {
    if(dynamicsD[k].id==id) {
      return(k);
      }
    }
   
  return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_identify_dynamic(SOLVER_t< complex<double> > & solver, int id, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  for(int k=0;k<dynamicsZ.size();k++) {
    if(dynamicsZ[k].id==id) {
      return(k);
      }
    }
   
  return(-1);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_initialize_dynamic(SOLVER_t<double> & solver, int id, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  for(int k=0;k<dynamicsD.size();k++) {
    if(dynamicsD[k].id==id) {
      solver.factorize=dynamicsD[k].factorize;
      solver.solve=dynamicsD[k].solve;
      solver.initialize=dynamicsD[k].initialize;
      solver.terminate=dynamicsD[k].terminate;
      solver.factorize=dynamicsD[k].factorize;
      solver.error=dynamicsD[k].error;
      solver.id=dynamicsD[k].id;
      }
    }
   
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_initialize_dynamic(SOLVER_t< complex<double> > & solver, int id, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  for(int k=0;k<dynamicsZ.size();k++) {
    if(dynamicsZ[k].id==id) {
      solver.factorize=dynamicsZ[k].factorize;
      solver.solve=dynamicsZ[k].solve;
      solver.initialize=dynamicsZ[k].initialize;
      solver.terminate=dynamicsZ[k].terminate;
      solver.factorize=dynamicsZ[k].factorize;
      solver.error=dynamicsZ[k].error;
      solver.id=dynamicsZ[k].id;
      }
    }
   
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int solver_register_dynamic_template(const char *name, vector<SOLVER_REG_t<T> > & dynamics,
                                                   int (*factorize)  (SOLVER_t<T> *, int, bool),
                                                   int (*solve)      (SOLVER_t<T> *, T *, T *, int, int, bool),
                                                   int (*initialize) (SOLVER_t<T> &, int, bool),
                                                   int (*terminate)  (SOLVER_t<T> *, int, bool),
                                                   int (*error)      (int, int, bool),
                                                   int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int id;
//   bool debug=false;
  
/*------------------------------------------------------------------------------
  temporary, until dynamic solver handling fully implemented */
//   for(int k=0;k<dynamics.size();k++) {
//     if(strcmp(name, dynamics[k].name)==0) {
//       TRAP_ERR_EXIT(-1, "solver %s already registered", name);
//       }
//     }
  
  SOLVER_REG_t<T> *reg=new SOLVER_REG_t<T>;
  
  reg->name=strdup(name);
  reg->factorize=factorize;
  reg->solve=solve;
  reg->initialize=initialize;
  reg->terminate=terminate;
  reg->factorize=factorize;
  reg->error=error;
  
/*------------------------------------------------------------------------------
  temporary, until dynamic solver handling fully implemented */
//   reg->id=100+dynamics.size();
  reg->id=get_SolverId(name, false);
  if(reg->id==-1) TRAP_ERR_EXIT(-1, "solver %s not known", name);
  
  dynamics.push_back(*reg);
  
  return(reg->id);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_register_dynamic_d(const char *name,
                                int (*factorize)  (SOLVER_t<double> *, int, bool),
                                int (*solve)      (SOLVER_t<double> *, double *, double *, int, int, bool),
                                int (*initialize) (SOLVER_t<double> &, int, bool),
                                int (*terminate)  (SOLVER_t<double> *, int, bool),
                                int (*error)      (int, int, bool),
                                int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int id;
  
  id=solver_register_dynamic_template(name, dynamicsD, factorize, solve, initialize, terminate, error, verbose, debug);
  
  return(id);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_register_dynamic_z(const char *name,
                                int (*factorize)  (SOLVER_t< complex<double> > *, int, bool),
                                int (*solve)      (SOLVER_t< complex<double> > *, complex<double> *, complex<double> *, int, int, bool),
                                int (*initialize) (SOLVER_t< complex<double> > &, int, bool),
                                int (*terminate)  (SOLVER_t< complex<double> > *, int, bool),
                                int (*error)      (int, int, bool),
                                int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int id;
  
  id=solver_register_dynamic_template(name, dynamicsZ, factorize, solve, initialize, terminate, error, verbose, debug);
  
  return(id);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solver_entry_dynamic(const char *name, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int id;
  
  SOLVER_REG_t<double>           tmpD;
  SOLVER_REG_t<complex<double> > tmpZ;
  
  id=solver_register_dynamic_template(name, dynamicsD, tmpD.factorize, tmpD.solve, tmpD.initialize, tmpD.terminate, tmpD.error, verbose, debug);
  id=solver_register_dynamic_template(name, dynamicsZ, tmpZ.factorize, tmpZ.solve, tmpZ.initialize, tmpZ.terminate, tmpZ.error, verbose, debug);
   
  return(id);
}





