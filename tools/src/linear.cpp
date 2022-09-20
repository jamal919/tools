
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

*******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <time.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "fe.h"
#include "matrix.h"

extern int printOK;
extern double ITER;
// int nbsyslinear=0;
const double diagonal=1.0;

extern char *get_solver_name(int solver_id);

/*----------------------------------------------------------------------------*/
/// gives a copy of an ordering
/**
by allocating independent arrays
whose values are initialised with that of the original.
Only the shared property is initialised to 0.
\date created 30 May 2011
\author Damien Allain

\param s source
\return an independent ordering
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

ordering_t *copy(const ordering_t *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ordering_t *r=new ordering_t;//returned value
  r->incidence=copy(s->incidence,s->pointer[s->nrows]);
  r->cardinal=copy(s->cardinal,s->nrows);
  r->pointer=copy(s->pointer,s->nrows+1);
  r->shared=0;
  r->nrows=s->nrows;
  r->max_cardinal=s->max_cardinal;
  r->hbw=s->hbw;
  r->type=s->type;
  return r;
}

/*----------------------------------------------------------------------------*/
/// gives a copy of a matrix
/**
by allocating independent packed and ordering.
\date created 30 May 2011
\author Damien Allain
\todo a copy_matrix_template for hyperzmatrix_t and hypermatrix_t is possible

\param s source
\return an independent matrix
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  hyperzmatrix_t copy(hyperzmatrix_t s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  hyperzmatrix_t r;//returned value
  if(s.neq()==-1)
    check_error(-1, "badly initialised source", __LINE__, __FILE__, 1);
  r.packed=copy(s.packed,s.ordering->pointer[s.neq()]);
  r.ordering=copy(s.ordering);
  r.shared_ordering=0;//means it is not shared
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class M_t,typename T> int sparsePPM_template_(M_t A,char *path,int d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *f;
  int i,j;//row and column indexes
  int n,ncols=-1;//value index
  T v;//value
  ordering_t *o=A.ordering;//shortcut
  if(o==NULL){
    fprintf(stderr,"giving up saving %s : null ordering\n",path);
    return 1;
    }
  f=fopen(path,"w");
  if(f==NULL){
    fprintf(stderr,"failed to open %s, ",path);
    path=strrchr(path,'/');path++;//get base name
    fprintf(stderr,"trying to save to %s instead\n",path);
    f=fopen(path,"w");
    }
  if(f==NULL){
    fprintf(stderr,"failed to open %s\n",path);
    return 1;
    }
  for(i=0;i<A.neq();i++){
    n=o->pointer[i+1]-1;
    if(ncols<o->incidence[n]){
      ncols=o->incidence[n];
      }
    }
  ncols++;
  if(d<1){
    d=1;
    }
  fprintf(f,"P6\n%d %d\n255\n",(int)ceil((double)ncols/d),(int)ceil((double)A.neq()/d));
  for(i=0;i<A.neq();i+=d){
    n=o->pointer[i];
    for(j=0;j<ncols;j+=d){
      for(;o->incidence[n]<j && n<o->pointer[i+1];n++);
      v=o->incidence[n]==j?A.packed[n]:0.;
      putc(log2(abs(v))*4+128,f);
      putc(log2(abs(v))*4+128,f);
      putc(log2(abs(v))*4+128,f);
      }
    }
  fclose(f);
}

/*----------------------------------------------------------------------------*/
///for use ONLY by sparsePPM()
/**
Do NOT call directly. Use sparsePPM() instead.
\date created 9 Jun 2011
\author Damien Allain
\todo color image
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int sparsePPM_(hypermatrix_t A,char *path,int d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=sparsePPM_template_<hypermatrix_t,double>(A,path,d);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int sparsePPM_(hyperzmatrix_t A,char *path,int d)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=sparsePPM_template_<hyperzmatrix_t,dcomplex>(A,path,d);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int LinearSystem_InitTriplet_template(T *M, bool duplicate, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,n,status;

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  26/02/2011:

    Triplet initialisation have been deported to calling routines, not
    fully checked

  KILLER:
    At the moment, real-valued system still relies on false CSC ordering,
    compensated by transposed system solving (hidden in poc-solver)

  KILLER:
    MUMPS modifies M->t.i and M->t.j
    
  KILLER:
    HIPS modifies M->t.i, M->t.j and M->packed
    
  KILLER:
    temporary patch will create memory leaks
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

/**----------------------------------------------------------------------------
  init triplet*/

  M->t.comptype=M->ordering->type;
  M->t.base=0;

  M->t.stype=0;

  if(M->neq()==0) {
    check_error(-1, "neq field is zero", __LINE__, __FILE__, 1);
    }
  M->t.nnz=M->ordering->pointer[M->neq()];
  M->t.nrow=M->neq();
  M->t.ncol=M->neq();

  if(duplicate) {
/**----------------------------------------------------------------------------
    mainly usefull in MPI mode where index are modified and CCO ordering used*/
    M->t.i=(int *) malloc( (M->neq()+1)*sizeof(int));
    M->t.j=(int *) malloc( M->ordering->pointer[M->neq()]*sizeof(int));
    for(n=0;n<M->neq()+1;n++) {
      M->t.i[n]=M->ordering->pointer[n];
      }
    for(m=0;m<M->ordering->pointer[M->neq()];m++) {
      M->t.j[m]=M->ordering->incidence[m];
      }
    }
  else {
    M->t.i = M->ordering->pointer;
    M->t.j = M->ordering->incidence;
    }
  
  M->t.x = M->packed;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_InitTriplet(hypermatrix_t *M, bool duplicate, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,n,status;
  status=LinearSystem_InitTriplet_template(M, duplicate, verbose, debug);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_InitTriplet(hyperzmatrix_t *M, bool duplicate, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int m,n,status;
  status=LinearSystem_InitTriplet_template(M, duplicate, verbose, debug);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_AdjustTriplet(hypermatrix_t *M, int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

/**----------------------------------------------------------------------------
  re-arrange matrix following solver type*/
  switch (solver_id) {
/*------------------------------------------------------------------------------
    DOMESTIC*/
    case SOLVER_ID_DOMESTIC:
      status=0;
      break;

/*------------------------------------------------------------------------------
    LAPACKF_*/
    case SOLVER_ID_LAPACK:
      status=0;
      break;

/*------------------------------------------------------------------------------
    SUNPERF*/
    case SOLVER_ID_SUNPERF:
      status=0;
      break;

/*------------------------------------------------------------------------------
    UMFPACK*/
    case SOLVER_ID_UMFPACK:
      status=0;
      break;
 
/*------------------------------------------------------------------------------
    MUMPS ASSYMETRIC*/
    case SOLVER_ID_MUMPS:
//      status=matrix_CSR2CSC(M);
      break;
 
/*------------------------------------------------------------------------------
    MUMPS SYMETRIC*/
    case SOLVER_ID_MUMPS_SYM:
      status=0;
      break;

/*------------------------------------------------------------------------------
    HIPS */
    case SOLVER_ID_HIPS:
      break;

/*------------------------------------------------------------------------------
    HYPRE */
    case SOLVER_ID_HYPRE:
      break;

/*------------------------------------------------------------------------------
    SpDOMESTIC*/
    case SOLVER_ID_SpDOMESTIC:
      status=0;
      break;

/*------------------------------------------------------------------------------
    PASTIX*/
    case SOLVER_ID_PASTIX:
      status=0;
      break;

    default:
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_factorisation_sequential_obsolete(hypermatrix_t *M, int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, lda, status, ierr=0;
  int *colstr, *rowind;
  char mtxtyp[2], cpivot, ordmthd;
  int outunt, msglvl = 0, error, nrhs, ldrhs, neqns, nz;
  float *rhs;
  double handle[150];
  int hbw, neq;

  char c;
  int *Tj;
  FILE *out;
  clock_t d1, d2;
  double duree1, duree2;
  int count, m, k, n;
  char **args;
  triplet *triplet;
  solver_t *solveur1;
  int verbose=0;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  Development notes

  Check :

  Notes 29/04/2007

    stype=0: unsymetric
    stype>0: matrix is square and symmetric. Entries in the lower triangular
    stype<0: matrix is square and symmetric. En

  26/02/2011:

    Triplet initialisation have been deported to calling routines, not checked

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=LinearSystem_AdjustTriplet(M, solver_id);
  triplet = &(M->t);

  d1 = clock();
  
  switch (solver_id) {
/*------------------------------------------------------------------------------
    LAPACKF_*/
    case SOLVER_ID_LAPACK:
      hbw   = M->hbw;
      if(hbw==-1) {
        M->ordering->gethbw();
        hbw   = M->ordering->hbw;
        }
      if(hbw==-1) {
        check_error(-1, "hbw (half-bandwidth) field is zero, check matrix construction routine...", __LINE__, __FILE__, 1);
        }
      printf("lapack band matrix solver\n");
#   ifdef LAPACKF_
      lda = 3 * hbw + 1;
      M->pivot = new int[M->neq()];
      status = band_unpack01(M->packed, M->ordering->incidence, M->ordering->cardinal, M->neq, hbw,&(M->band));
      dgbtrf_(&M->neq, &M->neq, &hbw, &hbw, M->band, &lda, M->pivot, &status);
      if(status!=0) {
        check_error(-1, "dgbtrf (double band matrix) factorization failed", __LINE__, __FILE__, 1);
        }
#   else
      check_error(-1, "please compile with -DLAPACKF_", __LINE__, __FILE__, 1);
#   endif
      break;
 
/*------------------------------------------------------------------------------
    SUNPERF*/
    case SOLVER_ID_SUNPERF:
      if(verbose==1) printf("sunperf packed matrix solver\n");
      rowind = M->ordering->incidence;
      colstr = M->ordering->pointer;
      neqns = neq;
      outunt = neq + 1;
      cpivot = 'n';
      strcpy(mtxtyp, "ss");
/*
      dgssin_ (mtxtyp,cpivot, &neqns, colstr, rowind, &outunt, &msglvl,
                 handle, &error);
      dgssor_ (ordmthd, handle, error);
      dgssfa_ (neqns, colstr, rowind, outunt, packed,handle, error);
      dgsssl_ (nrhs,rhs, ldrhs, handle, error);
*/
      break;

/*------------------------------------------------------------------------------
    UMFPACK*/
    case SOLVER_ID_UMFPACK:
      if(verbose==1) printf("UMFPACK matrix solver\n");
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("UMFPACK", triplet->stype);
      status = factorize(solveur1,triplet,verbose);
      if (status != 0) {
        cout << "UMFPACK factorize error:" << status << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      break;
 
/*------------------------------------------------------------------------------
    MUMPS ASSYMETRIC*/
    case SOLVER_ID_MUMPS:
#ifdef MUMPS
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("MUMPS", triplet->stype);
      ierr = factorize(solveur1,triplet,verbose);
      printf(" LinearSystem_factorisation (mono) status: %d\n",ierr);
      if (ierr != 0) {
        cout << "MUMPS factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
//      printf(" Symetrie = %d \n",coo_issym(M1));
      M->s  = solveur1;
#   else
      check_error(-1, "please compile with -DMUMPS", __LINE__, __FILE__, 1);
#   endif
      status=ierr;
      break;
 
/*------------------------------------------------------------------------------
    MUMPS SYMETRIC*/
    case SOLVER_ID_MUMPS_SYM:
#ifdef MUMPS
      ierr=0;
      triplet->stype=1; /*assume symetrie*/
      solveur1=(solver_t *) init_solver_obsolete("MUMPS", triplet->stype);
      ierr = factorize(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "MUMPS factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
//      printf(" Symetrie = %d \n",coo_issym(M1));
      M->s = solveur1;
#   else
      check_error(-1, "please compile with -DMUMPS", __LINE__, __FILE__, 1);
#   endif
      status=ierr;
      break;

/*------------------------------------------------------------------------------
    SpDOMESTIC*/
    case SOLVER_ID_SpDOMESTIC:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("SpDOMESTIC", triplet->stype);
      ierr = factorize(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "SpDOMESTIC factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      status=ierr;
      break;

/*------------------------------------------------------------------------------
    HIPS*/
    case SOLVER_ID_HIPS:
      if(verbose==1) printf("PASTIX matrix solver\n");
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("HIPS", triplet->stype);
      status = factorize(solveur1,triplet,verbose);
      if (status != 0) {
        cout << "HIPS factorize error:" << status << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      break;

/*------------------------------------------------------------------------------
    PASTIX*/
    case SOLVER_ID_PASTIX:
      if(verbose==1) printf("PASTIX matrix solver\n");
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("PASTIX", triplet->stype);
      status = factorize(solveur1,triplet,verbose);
      if (status != 0) {
        cout << "PASTIX factorize error:" << status << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      break;

    default:
      cout << "Solver_id = " << solver_id <<endl;
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;
    }
    
/**----------------------------------------------------------------------------
  store solver id */

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check : MANDATORY

  Notes

    If the following lines trigger a compiling fault, you must either:
    
      1) upgrade poc-solvers library to version 2 and successors and
         update T-UGOm Makefile consequently
      
      2) keep your usual poc-solvers version and add the "id" field
         in solver structure define in solverlib.h

         typedef struct solver_struct
         {
           char *name;
           int id;               <---------------------- here
           void *parameters;
           void *Mat;
         } solver;

         Re-compile library after change was made, the solverlib.h will
         be copied into the include directory.
         You might need to comment "#include mip.h" line in the
         include/solverlib.h as it is needed for library test program but
         redundant (and incomplete) for T-UGOm
    
    1st option is recommended, option 2 to be used as the last bullet
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  if(M->s==0) M->s=new solver_t;
  M->s->id=solver_id;
  
  d2 = clock();
  duree1=d2-d1;

  if(verbose==1) printf("time spent for factorisation: %.1lf s\n",duree1/CLOCKS_PER_SEC);
 
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_factorisation_sequential(hypermatrix_t *M, int solver_id, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, lda, status, ierr=0;
  int *colstr, *rowind;
  char mtxtyp[2], cpivot, ordmthd;
  int outunt, msglvl = 0, error, nrhs, ldrhs, neqns, nz;
  float *rhs;
  double handle[150];
  int hbw, neq;

  char c;
  int *Tj;
  FILE *out;
  clock_t d1, d2;
  double duree1, duree2;
  int count, m, k, n;

  triplet *triplet;
  solver_t *solveur1;
  
  char *solver_name=get_solver_name(solver_id);

  status=LinearSystem_AdjustTriplet(M, solver_id);
  triplet = &(M->t);

  d1 = clock();
  
  if(solver_id==SOLVER_ID_MUMPS_SYM) {
/*-----------------------------------------------------------------------------
    stype=0: unsymetric
    stype>0: matrix is square and symmetric. Entries in the lower triangular
    stype<0: matrix is square and symmetric. Entries in the upper triangular
-----------------------------------------------------------------------------*/
    triplet->stype=1; /*assume symetrie*/
//     printf(" Symetrie = %d \n",coo_issym(M1));
    }

  if(solver_id==SOLVER_ID_MUMPS) {
    triplet->stype=0; /*assume non-symetrie*/
    }
    
  switch (solver_id) {
/*----------------------------------------------------------------------
    DOMESTIC, not handled by poc-solver interface */
//     case SOLVER_ID_DOMESTIC:
//       if(verbose==1) printf("DOMESTIC band-matrix solver (now obsolete)\n");
//       hbw=M->hbw;
//       if(hbw==-1) {
//         M->ordering->gethbw();
//         hbw=M->hbw=M->ordering->hbw;
//         }
//       if(hbw==-1) {
//         check_error(-1, "hbw (half-bandwidth) field is zero, check matrix construction routine...", __LINE__, __FILE__, 1);
//         }
//       status = band_unpack02(M->packed, M->ordering->incidence, M->ordering->cardinal, M->neq, M->hbw,&(M->band));
//       status = d_LUDV(M->band, M->neq, M->hbw);
//       break;

/*----------------------------------------------------------------------
    LAPACKF_, not handled by poc-solver interface*/
    case SOLVER_ID_LAPACK:
      hbw   = M->hbw;
      if(hbw==-1) {
        M->ordering->gethbw();
        hbw=M->hbw=M->ordering->hbw;
        }
      if(hbw==-1) {
        check_error(-1, "hbw (half-bandwidth) field is zero, check matrix construction routine...", __LINE__, __FILE__, 1);
        }
      printf("lapack band matrix solver\n");
#   ifdef LAPACKF_
      lda = 3 * hbw + 1;
      M->pivot = new int[M->neq];
      status = band_unpack01(M->packed, M->ordering->incidence, M->ordering->cardinal, M->neq, hbw,&(M->band));
      dgbtrf_(&M->neq, &M->neq, &hbw, &hbw, M->band, &lda, M->pivot, &status);
      if(status!=0) {
        check_error(-1, "dgbtrf (double band matrix) factorization failed", __LINE__, __FILE__, 1);
        }
#   else
      check_error(-1, "please compile with -DLAPACKF_", __LINE__, __FILE__, 1);
#   endif
      break;
 
/*----------------------------------------------------------------------
    SUNPERF, deprecated */
    case SOLVER_ID_SUNPERF:
      check_error(-1, "SUNPERF solver deprecated", __LINE__, __FILE__, 1);
//       if(verbose==1) printf("sunperf packed matrix solver\n");
//       rowind = M->ordering->incidence;
//       colstr = M->ordering->pointer;
//       neqns = neq;
//       outunt = neq + 1;
//       cpivot = 'n';
//       strcpy(mtxtyp, "ss");
//       dgssin_ (mtxtyp,cpivot, &neqns, colstr, rowind, &outunt, &msglvl, handle, &error);
//       dgssor_ (ordmthd, handle, error);
//       dgssfa_ (neqns, colstr, rowind, outunt, packed,handle, error);
//       dgsssl_ (nrhs,rhs, ldrhs, handle, error); 
      break;
 
    case SOLVER_ID_UMFPACK:
    case SOLVER_ID_MUMPS:
    case SOLVER_ID_MUMPS_SYM:
    case SOLVER_ID_SpDOMESTIC: 
    case SOLVER_ID_HIPS:
    case SOLVER_ID_PASTIX:
//     case SOLVER_ID_MAPHYS:
/*----------------------------------------------------------------------
      handled by poc-solver interface*/
      if(verbose==1) printf("%s matrix solver\n",solver_name);
      solveur1=(solver_t *) init_solver_obsolete(solver_name, triplet->stype, verbose);
      status = factorize(solveur1,triplet,verbose);
      if (status != 0) {
        cout << solver_name << " factorize error:" << status << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      break;

    default:
      cout << "solver_id = " << solver_id <<endl;
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;
    } 
    
/*------------------------------------------------------------------------------
  store solver id */
  if(M->s==0) M->s=(solver_t *) malloc(sizeof(solver_t));/* You MUST use malloc() because the poc-solvers are coded in C. */
  M->s->id=solver_id;
  
  d2 = clock();
  duree1=d2-d1;

  if(verbose==1) printf("time spent for factorisation: %.1lf s\n",duree1/CLOCKS_PER_SEC);
 
  free(solver_name);
  
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int LinearSystem_factorisation_parallel(hypermatrix_t *M, int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    stype=0: unsymetric
    stype>0: matrix is square and symmetric. Entries in the lower triangular
    stype<0: matrix is square and symmetric. En

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status, ierr=0;
  int verbose=0;
  int error;
  int neq;
  triplet *triplet;
  solver_t *solveur1;
  

  if(M->distributor==0) {
    check_error(-1, "node distribution structure not initialized", __LINE__, __FILE__, 1);
    }

//   status = MatrixLocal2Global(M, *(M->d));
  triplet = &(M->t);
  
  printf("LinearSystem_factorisation_parallel, solver id=%d\n",solver_id);
  
  switch (solver_id) {
/*------------------------------------------------------------------------------
    MUMPS ASSYMETRIC*/
    case SOLVER_ID_MUMPS:
#ifdef MUMPS
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("MUMPS", triplet->stype);
      ierr = factorize(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "MUMPS factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      printf(" Symetrie = %d \n",coo_issym(triplet));
      M->s  = solveur1;
#   else
      check_error(-1, "please compile with -DMUMPS", __LINE__, __FILE__, 1);
#   endif
      status=ierr;
      break;
  
/*------------------------------------------------------------------------------
    MUMPS SYMETRIC*/
    case SOLVER_ID_MUMPS_SYM:
#ifdef MUMPS
      ierr=0;
      triplet->stype=1; /*assume symetrie*/
      solveur1=(solver_t *) init_solver_obsolete("MUMPS", triplet->stype);
      ierr = factorize(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "MUMPS factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      printf(" Symetrie = %d \n",coo_issym(triplet));
      M->s = solveur1;
#   else
      check_error(-1, "please compile with -DMUMPS", __LINE__, __FILE__, 1);
#   endif
      status=ierr;
      break;
  
/*------------------------------------------------------------------------------
    SpDOMESTIC*/
    case SOLVER_ID_SpDOMESTIC:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("SpDOMESTIC", triplet->stype);
      ierr = factorize(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "SpDOMESTIC factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s = solveur1;
      status=ierr;
      break;
/*------------------------------------------------------------------------------
    HYPRE */
    case SOLVER_ID_HYPRE:
#ifdef HYPRE
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("HYPRE", triplet->stype);
      ierr = factorize(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "HYPRE factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      printf(" Symetrie = %d \n",coo_issym(triplet));
      M->s  = solveur1;
#   else
      check_error(-1, "please compile with -DHYPRE", __LINE__, __FILE__, 1);
#   endif
      status=ierr;
      break;
/*------------------------------------------------------------------------------
    HIPS */
    case SOLVER_ID_HIPS:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("HIPS", triplet->stype);
      ierr = factorize(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "HIPS factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s  = solveur1;
      status=ierr;
      break;
   
    case SOLVER_ID_PASTIX:
      ierr=0;
      solveur1=(solver_t *) init_solver_obsolete("PASTIX", triplet->stype);
      ierr = factorize(solveur1,triplet,verbose);
      if (ierr != 0) {
        cout << "PASTIX factorize error:" << ierr << endl;
        check_error(-1, "solver failure", __LINE__, __FILE__, 1);
        }
      M->s  = solveur1;
      status=ierr;
      break;
   
 /*--------------------------------------------------------------------*/
    default:
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
    }
   
/**----------------------------------------------------------------------------
  store solver id */
  if(M->s==0) M->s=new solver_t;
  M->s->id=solver_id;
  printf(" LinearSystem_factorisation (parallel) status: %d\n",ierr);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int LinearSystem_factorisation(hypermatrix_t *M, int solver_id, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, lda, status, ierr;
  if(M->distributor==0) {
#ifdef HAVE_MPI
    status = MPI_Barrier (  MPI_COMM_WORLD);
#endif
    check_error(-1, "node distribution structure not initialized", __LINE__, __FILE__, 1);
#ifdef HAVE_MPI
    status = MPI_Barrier (  MPI_COMM_WORLD);
#endif
    }
  
/**----------------------------------------------------------------------------
  special treatment for diagonal matrices */
  status=matrix_IsDiagonal(*M);
  if(M->IsDiagonal) {
    if(verbose==1) check_error(-1, "diagonal matrix issue, check recent changes", __LINE__, __FILE__, 0);
    for(size_t n=0; n< M->neq(); n++) {
      if(M->packed[n]!=0) {
/**----------------------------------------------------------------------------
        WARNING : 2013-06-25, inversion cancelled to avoid matrix modification*/
//        M->packed[n]=1./M->packed[n];
        }
      else {
        return(-1);
        }
      }
    return(0);
    }
    
  switch(M->distributor->context) {
    case SEQUENTIAL_COMPUTING:
      status=LinearSystem_factorisation_sequential(M,  solver_id, verbose, debug);
      break;
    default:
      printf("LinearSystem_factorization, context=%d\n",M->distributor->context);
      check_error(-1, "invalid solver context", __LINE__, __FILE__, 1);
      break;
    }
//  M->ordering->type=M->t.comptype;
  return (status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *get_solver_name(int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  char *name;

  switch (solver_id) {
/*----------------------------------------------------------------------
    DOMESTIC*/
    case SOLVER_ID_DOMESTIC:
      name=strdup(SOLVER_NAME_DOMESTIC);
      break;

/*----------------------------------------------------------------------
    LAPACKF_*/
    case SOLVER_ID_LAPACK:
      name=strdup(SOLVER_NAME_LAPACK);
      break;

/*----------------------------------------------------------------------
    SUNPERF*/
    case SOLVER_ID_SUNPERF:
      name=strdup(SOLVER_NAME_SUNPERF);
      break;

/*----------------------------------------------------------------------
    UMFPACK*/
    case SOLVER_ID_UMFPACK:
      name=strdup(SOLVER_NAME_UMFPACK);
      break;
 
/*----------------------------------------------------------------------
    MUMPS ASSYMETRIC*/
    case SOLVER_ID_MUMPS:
      name=strdup(SOLVER_NAME_MUMPS);
      break;
 
/*----------------------------------------------------------------------
    MUMPS SYMETRIC*/
    case SOLVER_ID_MUMPS_SYM:
      name=strdup(SOLVER_NAME_MUMPS_SYM);
      break;

/*----------------------------------------------------------------------
    HIPS */
    case SOLVER_ID_HIPS:
      name=strdup(SOLVER_NAME_HIPS);
      break;

/*----------------------------------------------------------------------
    HYPRE */
    case SOLVER_ID_HYPRE:
      name=strdup(SOLVER_NAME_HYPRE);
      break;

/*----------------------------------------------------------------------
    SpDOMESTIC*/
    case SOLVER_ID_SpDOMESTIC: 
      name=strdup(SOLVER_NAME_SpDOMESTIC);
      break;

/*----------------------------------------------------------------------
    PASTIX*/
    case SOLVER_ID_PASTIX:
      name=strdup(SOLVER_NAME_PASTIX);
      break;

/*----------------------------------------------------------------------
    MAPHYS*/
    case SOLVER_ID_MAPHYS:
      name=strdup(SOLVER_NAME_MAPHYS);
      break;

    default:
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;
    }
    
  return(name);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_solve_sequential(hypermatrix_t & M, double *B)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, i, ndiff, ierr, solver_id;
  int hbw,neqs;
  int job, nrhs = 1,bandwidth;
  char cjob = 'N';
  double *x, *tmp, maxdiff, maxdiff2, rnorm;
  void *Numeric;
  clock_t d1, d2;
  double duree1, duree2;
  triplet *M1;
  solver_t *solveur1;
  int *lign,*colo;
  double *val;
  char filerhs[30], filemat[30];
  int transposed;
  double *Bglob;
  char filename[50];
  FILE *out;
  bool showperf=0;

  //printf("Entree LinearSystem_solve \n");
#ifdef DEBUG_PARALLEL
#ifdef HAVE_MPI
      status = MPI_Barrier (MPI_COMM_WORLD);
#endif
#endif

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check : MANDATORY

  Notes

    If the following lines trigger a compiling fault, you must either:
    
      1) upgrade poc-solvers library to version 2 and successors and
         update T-UGOm Makefile consequently
      
      2) keep your usual poc-solvers version and add the "id" field
         in solver structure define in solverlib.h

         typedef struct solver_struct
         {
           char *name;
           int id;               <---------------------- here
           void *parameters;
           void *Mat;
         } solver;

         Re-compile library after change was made, the solverlib.h will
         be copied into the include directory.
         You might need to comment "#include mip.h" line in the
         include/solverlib.h as it is needed for library test program but
         redundant (and incomplete) for T-UGOm
    
    1st option is recommended, option 2 to be used as the last bullet
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

//  printf("------------ Solve Sequential   -----------------------------------\n");

  if(M.s==0)
    check_error(-1, "solver structure not properly initialized", __LINE__, __FILE__, 1);
  solver_id=M.s->id;
  
  M1 =  &(M.t);
  solveur1=M.s;
  hbw   = M.hbw;
  neqs  = M.t.nrow;

  if (printOK == 1) {
    for (i=0; i <neqs; i++) {
      printf("       LinearSystem_solver RHS[%d]=  %e\n",i,B[i]);
      }
    }

  switch (solver_id) {
    case SOLVER_ID_LAPACK:
      hbw   = M.hbw;
      if(hbw==-1) {
        M.ordering->gethbw();
        hbw   = M.ordering->hbw;
        }
      if(hbw==-1) {
        check_error(-1, "hbw (half-bandwidth) field is zero, check matrix construction routine...", __LINE__, __FILE__, 1);
        }
      neqs   = M.t.nrow;
      bandwidth = 3 * hbw + 1;
#ifdef LAPACKF_
      dgbtrs_(&cjob, &neqs, &hbw, &hbw, &nrhs, M.band, &bandwidth, M.pivot, B, &neqs, &status);
      if(status!=0) {
        check_error(-1, "dgbtrs (double band matrix) solver failed", __LINE__, __FILE__, 1);
        }
#endif
      break;

    case SOLVER_ID_SUNPERF:
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
      break;

    case SOLVER_ID_UMFPACK:
      switch(M.ordering->type) {
        case ORDER_CSR:
          transposed=1;
          break;
        case ORDER_CSC:
          transposed=0;
          break;
        default:
          check_error(-1, "invalid ordering type", __LINE__, __FILE__, 1);
        }
      d1 = clock();
      ierr = solve(solveur1,B,transposed,false);
      d2 = clock();
      duree1 = d2 - d1;
      if(showperf) printf("UMFPACK solving, elapsed time=%e\n",duree1/CLOCKS_PER_SEC);
      break;

    case SOLVER_ID_MUMPS:
      d1 = clock();
      ierr = solve(solveur1,B,0,false);
      d2 = clock();
      duree1 = (d2 - d1)/CLOCKS_PER_SEC;
      if(showperf) printf("MUMPS solver, ellapsed=%e\n",duree1);
      break;

    case SOLVER_ID_MUMPS_SYM:
      d1 = clock();
      ierr = solve(solveur1,B,0,false);
      d2 = clock();
      duree1 = d2 - d1;
      if(showperf) printf("MUMPS solver, ellapsed=%e\n",duree1);
      break;

    case SOLVER_ID_SpDOMESTIC:
      switch(M.ordering->type) {
        case ORDER_CSR:
          transposed=1;
//          check_error(-1, "unstable mode?", __LINE__, __FILE__, 0);
          break;
        case ORDER_CSC:
          transposed=0;
          break;
        default:
          check_error(-1, "invalid ordering type", __LINE__, __FILE__, 1);
        }
      d1 = clock();
      ierr = solve(solveur1,B,transposed,false);
      d2 = clock();
      duree1 = d2 - d1;
      if(showperf) printf("SpDOMESTIC solving, elapsed time=%e\n",duree1/CLOCKS_PER_SEC);
      break;

    case SOLVER_ID_PASTIX:
      switch(M.ordering->type) {
        case ORDER_CSR:
          check_error(-1, "PASTIX needs CSC matrices", __LINE__, __FILE__, 1);
          break;
        case ORDER_CSC:
          transposed=0;
          break;
        default:
          check_error(-1, "invalid ordering type", __LINE__, __FILE__, 1);
        }
      d1 = clock();
      ierr = solve(solveur1,B,1,false);
      d2 = clock();
      duree1 = d2 - d1;
      if(showperf) printf("PASTIX solving, elapsed time=%e\n",duree1/CLOCKS_PER_SEC);
      break;

    default:
      check_error(-1, "solver not implemented", __LINE__, __FILE__, 1);
    }

  if (printOK == 1) {
    for (i=0; i <neqs; i++) {
      ITER=ITER+1;
      printf("       %e LinearSystem_solver Sol[%d]=  %e\n",ITER,i,B[i]);
      }
    }
//   nbsyslinear++;
//   sprintf(filename,"RHS.%d.%d.txt",solver_id,nbsyslinear);
//   out=fopen(filename,"w");
//   for (i=0; i <neqs; i++)   fprintf(out,"      %d -> %e \n",i,B[i]);
//   fclose(out);
   

  return(ierr);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int LinearSystem_solve_parallel_obsolete(hypermatrix_t M, double *B, mesh_t *mesh)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status, i, ndiff, ierr, solver_id;
//   int hbw,neqs;
// 
//   int job, nrhs = 1,bandwidth;
//   char cjob = 'N';
//   double *x, *tmp, maxdiff, maxdiff2, rnorm;
// 
//   void *Numeric;
//   double d1, d2;
//   double duree1, duree2;
//   /*Vec biter,xiter; */
//   triplet *M1;
//   solver_t *solveur1;
// 
//   double *val;
// 
//   int transposed;
//   double *Bglob;
//   double *Bglobtmp;
//   int allocated=0;
// 
//   int nglob,nprocs,rank,bcl, nndes, idi, k, nloc, root;
//     
// /**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// 
//   Development notes
// 
//   Check : MANDATORY
// 
//   Notes
// 
//     If the following lines trigger a compiling fault, you must either:
//     
//       1) upgrade poc-solvers library to version 2 and successors and
//          update T-UGOm Makefile consequently
//       
//       2) keep your usual poc-solvers version and add the "id" field
//          in solver structure define in solverlib.h
// 
//          typedef struct solver_struct
//            {
//            char *name;
//            int id;               <---------------------- here
//            void *parameters;
//            void *Mat;
//            } solver;
// 
//          Re-compile library after change was made, the solverlib.h will
//          be copied into the include directory.
//          You might need to comment "#include mip.h" line in the
//          include/solverlib.h as it is needed for library test program but
//          redundant (and incomplete) for T-UGOm
//     
//     1st option is recommended, option 2 to be used as the last bullet
//     
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
// 
//   printf("------------ Solve Parallel   -----------------------------------\n");
// 
//   if(M.s==0)
//     check_error(-1, "solver structure not properly initialized", __LINE__, __FILE__, 1);
//   solver_id=M.s->id;
// 
//   M1 =  &(M.t);
//   solveur1=M.s;
//   printOK = 1;
//  
// 
//   switch (solver_id) {
//     /*.......................................................*/
//     /* ..............   MUMPS ................................*/
//     case SOLVER_ID_MUMPS:
// #ifdef HAVE_MPI
// /**----------------------------------------------------------------------------
//    Centralised RHS on host 0 */
//       nndes = mesh->nvtxs;
//       nloc=0;
//       /**----------------------------------------------------------------------------
// 	 nb elts localy resolved */  
//       for (i=0; i <nndes; i++) {
// 	if ( (mesh->vertices[i].rank) == 0) {
// 	  nloc++;
// 	}
//       }
//       /**----------------------------------------------------------------------------
// 	 Nombre d'elts total: sum of each domain nb */
//       status = MPI_Allreduce ( &nloc, &nglob, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD );
//       //printf("       LinearSystem_solve_par %d : Nb elt total=%d\n",rank,nglob);
//       Bglob    = new double[nglob];
//       Bglobtmp = new double[nglob];
//       allocated=1;
//       
//       /**----------------------------------------------------------------------------
// 	 initialise exchange vector from local values (solved nodes) */
//       //printf("        Vecteurs Globaux Version PARALLEL \n",nndes,nloc);
//       for (i=0; i <nglob; i++) {
// 	Bglob[i]=0.0;
// 	Bglobtmp[i]=0.0;
//       }
//       for (i=0; i <nndes; i++) {
// 	if ( (mesh->vertices[i].rank) == 0) {
// 	  idi = mesh->vertices[i].origin;
// 	  //printf(" ligne=%d   val=%e \n",idi,B[i]);
// 	  Bglobtmp[idi] = B[i];
// 	}
//       }
//       
//       /**----------------------------------------------------------------------------
// 	 global vector construction: ok if domain are disjoint (summation) */
//       status = MPI_Allreduce ( Bglobtmp, Bglob, nglob, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
// #else
//       Bglob = B;
//       allocated=0;
// #endif
// #ifdef HAVE_MPI
//       d1 = MPI_Wtime();
// #endif
//       if (printOK == 1) {
// 	for (i=0; i <nglob; i++) {
// 	  printf("       LinearSystem_solver RHS[%d]=  %e\n",i,Bglob[i]);
// 	}
//         }
// /**----------------------------------------------------------------------------
//      call solver_t */
//      ierr = solve(solveur1,Bglob,0);
// #ifdef HAVE_MPI
//       d2 = MPI_Wtime();
// #endif
//       duree1 = d2 - d1;
//       if (printOK == 1) printf("MUMPS parallel solver, ellapsed=%e\n",duree1);
// #ifdef HAVE_MPI
// /**----------------------------------------------------------------------------
//       broadcast distributed solution*/
//       root=0;
//       status=  MPI_Bcast( Bglob, nglob, MPI_DOUBLE, root, MPI_COMM_WORLD);
//       if (printOK == 1) {
// 	for (i=0; i <nglob; i++) {
//           ITER=ITER+1;
// 	  printf("       %e LinearSystem_solver Sol[%d]=  %e\n",ITER,i,Bglob[i]);
//           }
//         }
// /**----------------------------------------------------------------------------
//       update local solution vector (all noddes)*/
//       for (i=0; i <nndes; i++) {
//         idi = mesh->vertices[i].origin;
//         B[i] = Bglob[idi];
//         }
// #endif
//       break;
// 
//  /*.......................................................*/
//  /*. ...................... HIPS ..........................*/
//      case SOLVER_ID_HIPS:
//        ierr = solve(solveur1,B,0);
//        break;
//        
//    
//     default:
//       check_error(-1, "solver not implemented in parallel", __LINE__, __FILE__, 1);
//     }
//   if (allocated) {
//     delete[] Bglob;
//     delete[] Bglobtmp;
//     }
//   return(ierr);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int LinearSystem_solve_parallel(hypermatrix_t & M, double *B)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, i, m, n, ierr, solver_id;
  int neqs;
  double d1, d2;
  double duree1, duree2;
  triplet *M1;
  solver_t *solveur1;
  int transposed;
  double *Bglob,*Bglobtmp;
  int nglob,rank, nndes, nloc, root,bcl;
  int allocated=0;
  int gmin,gmax,gsize,proc_id;
  int HIPSRESONB;
  char filename[50];
  FILE *out;
  
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check : MANDATORY

  Notes

    If the following lines trigger a compiling fault, you must either:
    
      1) upgrade poc-solvers library to version 2 and successors and
         update T-UGOm Makefile consequently
      
      2) keep your usual poc-solvers version and add the "id" field
         in solver structure define in solverlib.h

         typedef struct solver_struct
           {
           char *name;
           int id;               <---------------------- here
           void *parameters;
           void *Mat;
           } solver;

         Re-compile library after change was made, the solverlib.h will
         be copied into the include directory.
         You might need to comment "#include mip.h" line in the
         include/solverlib.h as it is needed for library test program but
         redundant (and incomplete) for T-UGOm
    
    1st option is recommended, option 2 to be used as the last bullet
    
@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
#ifdef HAVE_MPI
   ierr = MPI_Comm_rank(MPI_COMM_WORLD, &proc_id);
#endif

  if(M.s==0)
    check_error(-1, "solver structure not properly initialized", __LINE__, __FILE__, 1);
  solver_id=M.s->id;

  M1 =  &(M.t);
  solveur1=M.s;
  nndes = M.distributor->ndof;
  nloc  = M.distributor->nsolved;
  nglob = M.distributor->ngdof;

  gmin  = M.distributor->gmin;
  gmax  = M.distributor->gmax;
  gsize = M.distributor->gsize;

//  out=fopen("RESONB","r");
//  fscanf(out,"%d",&HIPSRESONB);
//  fclose(out);
//  if (proc_id ==0 ) {
//   out=fopen("RESONB","w");
//   HIPSRESONB++;
//   fprintf(out,"%d\n",HIPSRESONB);
//   fclose(out);
//  }

  switch (solver_id) {
    
/**----------------------------------------------------------------------------
    MUMPS */
    case SOLVER_ID_MUMPS:
#ifdef HAVE_MPI
/**----------------------------------------------------------------------------
  Centralised RHS on host 0 */
  
//   status = MPI_Barrier(MPI_COMM_WORLD);
//   printf("cpu %d: nglob=%d nloc=%d nndes=%d\n",gCPU_ID,nglob,nloc,nndes);
//   status = MPI_Barrier(MPI_COMM_WORLD);

/**----------------------------------------------------------------------------
  initialise exchange vector from local values (solved nodes) */
//   Bglob    = new double[nglob];
//   Bglobtmp = new double[nglob];
//   allocated=1;
//   for (i=0; i <nglob; i++) {
//     Bglob[i]   =0.0;
//     Bglobtmp[i]=0.0;
//     }
  Bglob    = M.mpi_global;
  Bglobtmp = M.mpi_tmp;
  for (i=0; i <nglob; i++) {
    Bglob[i]   =0.0;
    }
  //sprintf(filename,"LienarSolve.MUMPS.solvedLIST.%d.txt",proc_id);
  //out=fopen(filename,"w");
  for (i=0; i <M.d->nsolved; i++) {
    n = M.d->solved[i];
//     m = M.d->gindex[n];
    m = M.d->gsolved[i];
    Bglobtmp[m] = B[n];
    //fprintf(out, "%d => %d, %lf\n",m,n, B[n]);
    }
  //fclose(out);
  //sprintf(filename,"LienarSolve.MUMPS.solved.1.%d.txt",proc_id);
  //out=fopen(filename,"w");
    for (i=0; i <M.d->nsolved; i++) {
      m = M.d->gsolved[i];
      //  fprintf(out, "%d => %d, %lf\n",i,m,Bglobtmp[m]);
    }
    //fclose(out);

    // sprintf(filename,"LienarSolve.MUMPS.RHSL.1.%d.txt",proc_id);
    //   out=fopen(filename,"w");
    //   for (bcl=0; bcl < nloc; bcl++) fprintf(out, "%d => %lf\n",bcl,B[bcl]);
    //   fclose(out);

/**----------------------------------------------------------------------------
  global vector construction: ok if domain are disjoint (summation) */
  status = MPI_Allreduce ( Bglobtmp, Bglob, nglob, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
//  status = MPI_Allreduce ( &(Bglobtmp[gmin]), &(Bglob[gmin]), gsize, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
#else
  allocated=0;
  Bglob = B;
#endif
#ifdef HAVE_MPI
      d1 = MPI_Wtime();
#endif
/**----------------------------------------------------------------------------
      call solver_t */
      ierr = solve(solveur1,Bglob,0,false);
#ifdef HAVE_MPI
      d2 = MPI_Wtime();
#endif
      duree1 = d2 - d1;
      printf("MUMPS parallel solver, ellapsed=%e\n",duree1);
#ifdef HAVE_MPI
/**----------------------------------------------------------------------------
      broadcast distributed solution*/
      root=0;
      status=  MPI_Bcast( Bglob, nglob, MPI_DOUBLE, root, MPI_COMM_WORLD);
/**----------------------------------------------------------------------------
      update local solution vector (all noddes)*/
      //sprintf(filename,"LienarSolve.MUMPS.Glob2locT.%d.txt",proc_id);
      //out=fopen(filename,"w");
      for (n=0; n <nndes; n++) {
        m = M.d->gindex[n];
        B[n] = Bglob[m];
        //fprintf(out, "%d => %d \n",m,n);
        }
      //fclose(out);
      
//      sprintf(filename,"MUMPS/LinearSolve.SOL.2.%d.%d.txt",HIPSRESONB,proc_id);
#endif
        break;

/*.......................................................*/
/*. ...................... HIPS ..........................*/
///  VERSION GLOBALE
  case SOLVER_ID_HIPS:
 #ifdef HAVE_MPI
    Bglob    = M.mpi_global;
    Bglobtmp = M.mpi_tmp;
    for (i=0; i <nglob; i++) Bglob[i]   =0.0;
   
    for (i=0; i <M.d->nsolved; i++) {
      n = M.d->solved[i];
      //     m = M.d->gindex[n];
      m = M.d->gsolved[i];
      Bglobtmp[m] = B[n];
    }
/**----------------------------------------------------------------------------
  global vector construction: ok if domain are disjoint (summation) */
  status = MPI_Allreduce ( Bglobtmp, Bglob, nglob, MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD );
  /// ......................... SOLVE ......
  ierr = solve(solveur1,Bglob,0); //.........................................
  
/**----------------------------------------------------------------------------
      update local solution vector (all noddes)*/
  //sprintf(filename,"LienarSolve.HIPS.Glob2locT.%d.txt",proc_id);
  //  out=fopen(filename,"w");
      for (n=0; n <nndes; n++) {
        m = M.d->gindex[n];
        B[n] = Bglob[m];
        //fprintf(out, "%d => %d \n",m,n);
        }
      //fclose(out);
//      sprintf(filename,"HIPS/LienarSolve.SOL.2.%d.%d.txt",HIPSRESONB,proc_id);
#endif
    break;
//  /*.......................................................*/
//  /*. ...................... HIPS ..........................*/
///  VERSION LOCALE
//      case SOLVER_ID_HIPS:
//        Bglob    = M.mpi_global;
//        Bglobtmp = M.mpi_tmp;
//        for (i=0; i <nglob; i++) {
// 	 Bglob[i]   =0.0;
//        }
//         sprintf(filename,"LienarSolve.HIPS.solvedLIST.%d.txt",proc_id);
//        out=fopen(filename,"w");
//        for (i=0; i <M.d->nsolved; i++) {
// 	 n = M.d->solved[i];
// 	 //     m = M.d->gindex[n];
// 	 m = M.d->gsolved[i];
// 	 Bglob[i] = B[n];
//  	 fprintf(out, "%d => %d, %lf\n",i,m,Bglob[m]);
//        }
//        fclose(out);

//        sprintf(filename,"LienarSolve.HIPS.RHSL.1.%d.txt",proc_id);
//        out=fopen(filename,"w");
//        for (bcl=0; bcl < nloc; bcl++) fprintf(out, "%d => %lf\n",bcl,Bglob[bcl]);
//        fclose(out);
//        /// ......................... SOLVE ......
//        ierr = solve(solveur1,Bglob,0); //.........................................
//        /// ......................... .... ......
//        sprintf(filename,"LienarSolve.HIPS.SOL.1.%d.txt",proc_id);
//        out=fopen(filename,"w");
//        for (n=0; n < nloc; n++) {
// 	 m = M.d->gindex[n];
// 	 B[n] = Bglob[n];
// 	 fprintf(out, "%d => %d, %lf\n",n,m,B[n]);
//        }
//        fclose(out);
//        sprintf(filename,"LienarSolve.HIPS.RES.solvedLIST.%d.txt",proc_id);
//        out=fopen(filename,"w");
//        for (i=0; i <M.d->nsolved; i++) {
// 	  n = M.d->solved[i];
// 	  //     m = M.d->gindex[n];
// 	  m = M.d->gsolved[i];
//           B[n]=Bglob[i];
//  	  fprintf(out, "%d => %d, %d, %lf\n",i,m,n,B[n]);
//        }
//         fclose(out);
// //         sprintf(filename,"LienarSolve.HIPS.RES.GLOB.%d.txt",proc_id);
// //         out=fopen(filename,"w");
// //         for (n=0; n <nndes; n++) {
// //            m = M.d->gindex[n];
// //            B[n] = Bglob[m];
// //            fprintf(out, "%d => %d, %d, %lf\n",i,m,n,B[n]);
// //         }
// //         fclose(out);

// //        out=fopen(filename,"w");
// //        for (bcl=0; bcl < nloc; bcl++) {
// // 	 m = M.d->gindex[bcl];
// // 	 fprintf(out, "%d => %d, %lf\n",bcl,m,B[bcl]);
// //        }
// //        fprintf(out, "MARQUE\n");
// //        fclose(out);
//         sprintf(filename,"LienarSolve.HIPS.SOL.2.%d.txt",proc_id);
       
//        break;

/**----------------------------------------------------------------------------
    OTHERS */
    default:
      printf("matrix %s\n",M.name.c_str());
      check_error(-1, "solver not implemented for parallel computing", __LINE__, __FILE__, 1);
      break;
    }
    
//   if (allocated) {
//     delete[] Bglob;
//     delete[] Bglobtmp;
//     }

//  out=fopen(filename,"w");
//  for (bcl=0; bcl < nloc; bcl++) fprintf(out, "%d => %lf\n",bcl,B[bcl]);
//   fclose(out);
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
//  TRAP_ERR_EXIT(-23,"exiting\n");
 
  return(ierr);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_solve(hypermatrix_t & M, double *B)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  extern int print_vector(char* file, mesh_t mesh, int discretisation, double *val);
  
  if(M.distributor==0) {
#ifdef HAVE_MPI
    status = MPI_Barrier(MPI_COMM_WORLD);
#endif
    check_error(-1, "node distribution structure not initialized", __LINE__, __FILE__, 1);
#ifdef HAVE_MPI
    status = MPI_Barrier(MPI_COMM_WORLD);
#endif
    }
    
//  printf("solve %s matrix with %s solver\n",M.name.c_str(),M.s->name);
  
/**----------------------------------------------------------------------------
  special treatment for diagonal matrices */
  if(M.IsDiagonal) {
    check_error(-1, "diagonal matrix issue, check recent changes", __LINE__, __FILE__, 0);
/**----------------------------------------------------------------------------
    WARNING : 2013-06-25, inversion cancelled to avoid matrix modification*/
    for(size_t n=0; n< M.neq(); n++) {
//      B[n]*=M.packed[n];
      B[n]/=M.packed[n];
      }
    return(0);
    }

  switch(M.distributor->context) {
    case SEQUENTIAL_COMPUTING:
//      printf("call LinearSystem_solve_sequential, context=%d\n",M.d->context);
      status= LinearSystem_solve_sequential(M, B);
      break;
    case PARALLEL_COMPUTING:
//      printf("call LinearSystem_solve_parallel_new, context=%d\n",M.d->context);
      status= LinearSystem_solve_parallel(M, B);
//      status= LinearSystem_solve_parallel(M, B, &(gFEmesh[0]));
      break;
    default:
      printf("LinearSystem_solve, context=%d\n",M.distributor->context);
      check_error(-1, "invalid solver context", __LINE__, __FILE__, 1);
      break;
    }
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_solve(hypermatrix_t & M, complex<double> *B)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,nMax=M.neq(),status;
  double*Br,*Bi;

  Br=new double[nMax];
  Bi=new double[nMax];

  for(n=0;n<nMax;n++){
    Br[n]=real(B[n]);
    Bi[n]=imag(B[n]);
    }
  status=LinearSystem_solve(M, Br);
  status=LinearSystem_solve(M, Bi);
  for(n=0;n<nMax;n++){
    B[n]=complex<double>(Br[n],Bi[n]);
    }
  delete[]Br;
  delete[]Bi;

  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int LinearSystem_print(char *file, hypermatrix_t A, double *b, double *x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int bcl;
  triplet t;
  FILE *fp;
  
/**----------------------------------------------------------------------
  Print System */
/**----------------------------------------------------------------------------
  MPI unsafe */
  fp = fopen(file, "w");
  t = A.t;
  fprintf(fp,"  %%%MatrixMarket matrix coordinate real general \n");
  fprintf(fp,"%d  %d  %d 1 0 \n",t.nrow,t.nrow,t.nnz);
  for (bcl=0; bcl < t.nnz; bcl++) {
    fprintf(fp,"%d  %d  %e \n", t.i[bcl],t.j[bcl],t.x[bcl]);
    }
  //fprintf(fp," Vecteur b\n");
  for (bcl=0; bcl < t.nrow; bcl++) {
    //fprintf(fp,"b[%d]=%e\n",bcl,b[bcl]);
    fprintf(fp,"%d  %e\n",bcl,b[bcl]);
    }
  return fclose(fp);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_initialise(hypermatrix_t *A, int *clamped, int nclamped, int solver_id, bool force_duplicate, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j, k, n, asymmetric=1;
  int reorder;
  ordering_t *CSCordering;
  int transposed;
  bool duplicate=false;
  int MPIsize;
  bool spokesman=(gCPU_ID==gCPU_MASTER);
  
  if(spokesman and verbose) printf("%s : cpu=%3d, matrix is < %s >\n",__FUNCTION__, gCPU_ID, A->name.c_str());
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  set the matrix row to 0---0 1 0---0 at listed nodes
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(k = 0; k < nclamped; k++) {
    n=clamped[k];
    for(j = A->ordering->pointer[n]; j < A->ordering->pointer[n + 1]; j++) {
      if(A->ordering->incidence[j] == n) {
        A->packed[j] = diagonal;
        } 
      else {
        A->packed[j] = 0.0;
        }
      }
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  periodic mesh implementation
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   status=PackedMatrix_periodic(gFEmesh[0], gPeriodic[0], *A, A->discretisation);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  Deal with matrix ordering, as expected by the external sparse matrix solvers
  
  PASTIX special case:
  
    PASTIX expects CSC ordering. This is an issue when doing distributed
    resolution (T-UGOm natural workout is to provide distributed CSR matrices).
    
    Using the property that a matrix with CSR ordering is identical to its 
    transpose matrix with the same ordering BUT flagged as being CSC, we provide
    then to PASTIX the native matrix, flagging it CSR ordering as being CSC and 
    requesting then the solution of the transpose system.
 
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#define NEW_PASTIX_WORKOUT 1

  switch (solver_id) {
    case SOLVER_ID_DOMESTIC:
    case SOLVER_ID_LAPACK:
    case SOLVER_ID_SUNPERF:
/**----------------------------------------------------------------------
      not handled by poc-solvers, do nothing */
      break;
    case SOLVER_ID_HIPS:
    case SOLVER_ID_MUMPS:
/**----------------------------------------------------------------------
      leave HIPS in CSR numbering, necessary to build a proper COO (row-arranged) */
      asymmetric=0;  // HERE!!!
      break;
    case SOLVER_ID_PASTIX:
#if NEW_PASTIX_WORKOUT
/*------------------------------------------------------------------------------
      passing transpose matrix*/
      A->ordering->type=ORDER_CSC;
      transposed=1;
      break;
#endif
    case SOLVER_ID_UMFPACK:
    case SOLVER_ID_SpDOMESTIC:
/**----------------------------------------------------------------------
      convert to CSC numbering, necessary to build a proper COO (column-arranged) */
      asymmetric= matrix_check_symmetry(*A->ordering);
      if(asymmetric==0) {
/**----------------------------------------------------------------------
        if ordering is symmetric, CSR and CSC ordering are the same
        thus limit operation to matrix geometrical transposition */
        reorder=0;
        status= matrix_CSR2CSC(A);
        }
      else { 
/**----------------------------------------------------------------------
        if ordering is not symmetric, we need a separate CSC ordering
        in addition to matrix geometrical transposition */
        reorder=1;
        status=matrix_check_Zero(*A);
        CSCordering=ordering_CSR2CSC(A->ordering);
        status= matrix_CSR2CSC(A,CSCordering);
        A->ordering->destroy();
        A->ordering=CSCordering;
        status= matrix_check_consistency(*CSCordering);
        status=matrix_check_Zero(*A);
        }
      break; 
    default:
      TRAP_ERR_EXIT(-1,"%s cpu=%d : unknown solver\n", __func__, gCPU_ID);
      break;
    } 

  if(A->distributor==0) {
    check_error(-1, "node distribution structure not initialized", __LINE__, __FILE__, 1);
    }

#ifdef HAVE_MPI
  MPIsize = MPI::COMM_WORLD.Get_size();
#else
  MPIsize = 1;
#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  Solvers might modify or modify AND re-allocate structure arrays (pointer and incidence)
  as well as matrix itself (packed array).
  
  There are 3 consequences:
  
  * modified ordering not yet understood/handled in matrix routines
  
  * matrix values modified (hence troubles if re-using)
  
  * memory de-allocation messed up

  In case of non MPI it should be possible to manage the full problem without duplicating
  arrays.

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  initialise triplet, temporary patch for PARALLEL runs 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  switch(A->distributor->context) {
    case SEQUENTIAL_COMPUTING:
/**----------------------------------------------------------------------
      in case of non-symmetric matrices, ordering and matrix will be modified by solver */
      if(force_duplicate) duplicate=true;
      else duplicate=(asymmetric!=0);
      if(verbose==1 and spokesman) printf("%s: matrix %s, system duplication=%d\n", __FUNCTION__, A->name.c_str(), duplicate);
      status=LinearSystem_InitTriplet(A,duplicate,verbose,debug);
      break;
    case PARALLEL_COMPUTING:
      duplicate=true;
      if(verbose==1 and spokesman) printf("%s: matrix %s, system duplication=%d\n", __FUNCTION__, A->name.c_str(), duplicate);
      status=LinearSystem_InitTriplet(A,duplicate,verbose,debug);
      break;
    default:
      printf("LinearSystem_factorization, context=%d\n",A->distributor->context);
      check_error(-1, "invalid solver context", __LINE__, __FILE__, 1);
      break;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  factorize matrix 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status = LinearSystem_factorisation(A, solver_id, verbose, debug);
  
  if(A->s!=0) A->s->transpose_matrix=transposed;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  additional checks
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(A->packed!=A->t.x) {
    if(debug and spokesman) check_error(-1, "matrix ("+A->name+") has been re-allocated by solver", __LINE__, __FILE__, 0);
    if(duplicate==false) {
      A->packed=A->t.x; 
      if(debug) check_error(-1, "patch applies", __LINE__, __FILE__, 0);
      }
    }
    
  if(A->ordering->pointer!=A->t.i) {
    if(debug and spokesman) check_error(-1, "pointer arrray ("+A->name+") has been re-allocated by solver", __LINE__, __FILE__, 0);
    if(duplicate==false) {
      A->ordering->pointer=A->t.i;
      if(debug) check_error(-1, "patch applies", __LINE__, __FILE__, 0);
      }
    }
  
  if(A->ordering->incidence!=A->t.j) {
    if(debug and spokesman) check_error(-1, "incidence array ("+A->name+") has been re-allocated by solver", __LINE__, __FILE__, 0);
    if(duplicate==false) {
      A->ordering->incidence=A->t.j;
      if(debug) check_error(-1, "patch applies", __LINE__, __FILE__, 0);
      }
    }  
  
  if(A->ordering->type!=A->t.comptype) {
    if(debug and spokesman) check_error(-1, "ordering type has been changed by solver", __LINE__, __FILE__, 0);
    if(duplicate==false) {
      A->ordering->type=A->t.comptype;
      if(debug) check_error(-1, "patch applies", __LINE__, __FILE__, 0);
      }
    }  
  if(status!=0) printf("Real-valued linear system factorisation status: %d\n", status);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int LinearSystem_initialise_obsolete(hypermatrix_t *A, int *clamped, int nclamped, int solver_id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int j, k, n;
  int reorder;
  ordering_t *CSCordering;
  int verbose=0;
  bool debug=false;

/**----------------------------------------------------------------------
  set the matrix row to 0---0 1 0---0 at listed nodes*/
  for(k = 0; k < nclamped; k++) {
    n=clamped[k];
    for(j = A->ordering->pointer[n]; j < A->ordering->pointer[n + 1]; j++) {
      if(A->ordering->incidence[j] == n) {
        A->packed[j] = diagonal;
        }
      else {
        A->packed[j] = 0.0;
        }
      }
    }
    
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  24/10/2010:

    KILLER : CSC expect column sparse matrix, we provide row sparse matrix
             up to now, this was overcome by using the transpose solver, it
             does work as our matrices are geometrically symmetric

  01/02/2012:

    some proper matrix transformation efforts done. Interaction with
    poc-solvers to be monitored...
    At the moment, try to furnish CSC matrices, as it is the natural
    ordering for some solvers, and poc-solvers provides interfaces
    to convert from CSC to the other ordering.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  switch (solver_id) {
    case SOLVER_ID_DOMESTIC:
    case SOLVER_ID_LAPACK:
    case SOLVER_ID_SUNPERF:
//    case SOLVER_ID_HIPS:
/**----------------------------------------------------------------------
      HIPS needs COO structure, handled by poc-solvers*/
      break;
    default:
      status= matrix_check_symmetry(*A->ordering);
      if(status==0) {
/**----------------------------------------------------------------------
        if ordering is symmetric, CSR and CSC ordering are the same*/
        reorder=0;
        status= matrix_CSR2CSC(A);
        }
      else {
/**----------------------------------------------------------------------
        if ordering is not symmetric, wee need a separate CSC ordering*/
        reorder=1;
        status=matrix_check_Zero(*A);
        CSCordering=ordering_CSR2CSC(A->ordering);
        status= matrix_CSR2CSC(A,CSCordering);
        A->ordering->destroy();
        A->ordering=CSCordering;
        status= matrix_check_consistency(*CSCordering);
        status=matrix_check_Zero(*A);
        }
      break;
    }

/**----------------------------------------------------------------------
  initialise triplet, temporary patch for PARALLEL runs */
  switch(A->distributor->context) {
    case SEQUENTIAL_COMPUTING:
      status=LinearSystem_InitTriplet(A,false,0,false);
      break;
    case PARALLEL_COMPUTING:
/**----------------------------------------------------------------------
      duplicate arrays, else sparse matrix operations garbled */
      status=LinearSystem_InitTriplet(A,true,0,false);
      break;
    default:
      printf("LinearSystem_factorization, context=%d\n",A->distributor->context);
      check_error(-1, "invalid solver context", __LINE__, __FILE__, 1);
      break;
    }

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status = LinearSystem_factorisation(A, solver_id, verbose, debug);

//   if(A->ordering->incidence[0]!=0) {
//     printf("\nwarning, %s numbering changed by solver %d\n\n",A->name.c_str(),solver_id);
//     }

  if(status!=0) printf("Real-valued linear system factorisation status: %d\n", status);

  return (status);
}

