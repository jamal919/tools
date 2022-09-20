
/*******************************************************************************

  poc_solvers, SIROCCO 2006-2018

  Unstructured Ocean Grid initiative

Developers:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
Contributors:
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France

*******************************************************************************/

#include <string.h>
#include <omp.h>
#include "poc-solvers.h"
#include "solvers-interface.h"

int HIPS_SOLVERNB =0;
int HIPS_SOLVERNBZ=0;
  
static int *DHIPS_SOLVER=0, *ZHIPS_SOLVER=0;

#define HIPS_PARTITION_TYPE dHIPS_PARTITION_TYPE
#define HIPS_DOMSIZE dHIPS_DOMSIZE
#define HIPS_ITERATIVE dHIPS_ITERATIVE

#define HIPS_DROPTOL0   dHIPS_DROPTOL0
#define HIPS_DROPTOL1   dHIPS_DROPTOL1
#define HIPS_PREC       dHIPS_PREC

#define HIPS_VERBOSE dHIPS_VERBOSE
#define HIPS_SYMMETRIC dHIPS_SYMMETRIC
#define HIPS_FORTRAN_NUMBERING dHIPS_FORTRAN_NUMBERING
    
#define HIPS_ASSEMBLY_OVW dHIPS_ASSEMBLY_OVW
#define HIPS_ASSEMBLY_FOOL dHIPS_ASSEMBLY_FOOL
  
  
// #define  d
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  wrappers necessary to template use
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#ifdef HIPS

INTS HIPS_SetGlobalRHS(INTS id, double *b, INTS proc_root, INTS op) {

  int status;
  status=dHIPS_SetGlobalRHS(id, b, proc_root, op);
  return(status);
}

INTS HIPS_SetGlobalRHS(INTS id, complex<double> *b, INTS proc_root, INTS op) {

  int status;
  status=zHIPS_SetGlobalRHS(id, b, proc_root, op);
  return(status);
}


INTS HIPS_SetRHS(INTS id, INTS unknownnbr, INTS *unknownlist, double *b, INTS op, INTS op2, INTS mode)  {

  int status;
  status=dHIPS_SetRHS(id, unknownnbr, unknownlist, b, op, op2, mode);
  return(status);
}

INTS HIPS_SetRHS(INTS id, INTS unknownnbr, INTS *unknownlist, complex<double> *b, INTS op, INTS op2, INTS mode)  {

  int status;
  status=zHIPS_SetRHS(id, unknownnbr, unknownlist, b, op, op2, mode);
  return(status);
}

INTS HIPS_GetGlobalSolution(INTS id, double *x, INTS root) {

  int status;
  status=dHIPS_GetGlobalSolution(id,  x, root);
  return(status);
} 

INTS HIPS_GetGlobalSolution(INTS id, complex<double> *x, INTS root) {

  int status;
  status=zHIPS_GetGlobalSolution(id,  x, root);
  return(status);
} 

INTS HIPS_GetSolution(INTS id,  INTS n, INTS *nlist, double *x, INTS mode) {

  int status;
  status=dHIPS_GetSolution(id,  n, nlist, x, mode);
  return(status);
}

INTS HIPS_GetSolution(INTS id,  INTS n, INTS *nlist, complex<double> *x, INTS mode) {

  int status;
  status=zHIPS_GetSolution(id,  n, nlist, x, mode);
  return(status);
}

INTS HIPS_SetOptionINT(hips_t<double> *parameters, INTS id, INTS number, INTS value) {

  int status;
  status=dHIPS_SetOptionINT(id, number, value);
  return(status);
}
  
INTS HIPS_SetOptionINT(hips_t< complex<double> > *parameters, INTS id, INTS number, INTS value) {

  int status;
  status=zHIPS_SetOptionINT(id, number, value);
  return(status);
};


INTS HIPS_SetOptionREAL(hips_t<double> *parameters, INTS id, INTS number, REAL value) {

  int status;
  status=dHIPS_SetOptionREAL(id, number, value);
  return(status);
}
  
INTS HIPS_SetOptionREAL(hips_t< complex<double> > *parameters, INTS id, INTS number, REAL value) {

  int status;
  status=zHIPS_SetOptionREAL(id, number, value);
  return(status);
};


INTS HIPS_SetDefaultOptions(hips_t<double> *parameters, INTS id, INTS stratnum) {

  int status;
  status=dHIPS_SetDefaultOptions(id, stratnum);
  return(status);
}
  
INTS HIPS_SetDefaultOptions(hips_t< complex<double> > *parameters, INTS id, INTS stratnum) {

  int status;
  status=zHIPS_SetDefaultOptions(id, stratnum);
  return(status);
};




INTS HIPS_AssemblyBegin(SOLVER_t<double> *solver, INTS id, INTL nnz, INTS op, INTS op2, INTS mode, INTS symmetric) {

  int status;
  status=dHIPS_AssemblyBegin(id, nnz, op, op2, mode, symmetric);
  return(status);
};

INTS HIPS_AssemblyBegin(SOLVER_t< complex<double> > *solver, INTS id, INTL nnz, INTS op, INTS op2, INTS mode, INTS symmetric) {

  int status;
  status=zHIPS_AssemblyBegin(id, nnz, op, op2, mode, symmetric);
  return(status);
};

INTS HIPS_AssemblySetValue(INTS id, INTS row, INTS col, double value) {

  int status;
  status=dHIPS_AssemblySetValue(id, row, col, value);
  return(status);
};

INTS HIPS_AssemblySetValue(INTS id, INTS row, INTS col, complex<double> value) {

  int status;
  status=zHIPS_AssemblySetValue(id, row, col, value);
  return(status);
};

INTS HIPS_AssemblyEnd(SOLVER_t<double> *solver, INTS id) {

  int status;
  status=dHIPS_AssemblyEnd(id);
  return(status);
};

INTS HIPS_AssemblyEnd(SOLVER_t< complex<double> > *solver, INTS id) {

  int status;
  status=zHIPS_AssemblyEnd(id);
  return(status);
};


void HIPS_ExitOnError(hips_t<double> *parameters, INTS ierr) {

  int status;
  dHIPS_ExitOnError(ierr);
};

void HIPS_ExitOnError(hips_t< complex<double> > *parameters, INTS ierr) {

  int status;
  zHIPS_ExitOnError(ierr);
};
    
  
INTS HIPS_SetCommunicator(hips_t<double> *parameters, INTS id, MPI_Comm mpicom) {

  int status;
  status=dHIPS_SetCommunicator(id, mpicom);
  return(status);
};

INTS HIPS_SetCommunicator(hips_t< complex<double> > *parameters, INTS id, MPI_Comm mpicom) {

  int status;
  status=zHIPS_SetCommunicator(id, mpicom);
  return(status);  
};


INTS HIPS_SetPartition(hips_t<double> *parameters, INTS id, INTS ndom, INTS *mapptr, INTS *mapp) {

  int status;
  status=dHIPS_SetPartition(id, ndom, mapptr, mapp);
  return(status);
};

INTS HIPS_SetPartition(hips_t< complex<double> > *parameters, INTS id, INTS ndom, INTS *mapptr, INTS *mapp) {

  int status;
  status=dHIPS_SetPartition(id, ndom, mapptr, mapp);
  return(status);  
};


INTS HIPS_GraphBegin(hips_t<double> *parameters, INTS id, INTS n, INTL edgenbr) {

  int status;
  status=dHIPS_GraphBegin(id, n, edgenbr);
  return(status);
};

INTS HIPS_GraphBegin(hips_t< complex<double> > *parameters, INTS id, INTS n, INTL edgenbr) {

  int status;
  status=zHIPS_GraphBegin(id, n, edgenbr);
  return(status);  
};


INTS HIPS_GraphEdge(hips_t<double> *parameters, INTS id, INTS col, INTS row) {

  int status;
  status=dHIPS_GraphEdge(id, col, row);
  return(status);
};

INTS HIPS_GraphEdge(hips_t< complex<double> > *parameters, INTS id, INTS col, INTS row) {

  int status;
  status=zHIPS_GraphEdge(id, col, row);
  return(status);  
};


INTS HIPS_GraphEnd(hips_t<double> *parameters, INTS id) {

  int status;
  status=dHIPS_GraphEnd(id);
  return(status);
};

INTS HIPS_GraphEnd(hips_t< complex<double> > *parameters, INTS id) {

  int status;
  status=zHIPS_GraphEnd(id);
  return(status);  
};


INTS HIPS_Clean(hips_t<double> *parameters, INTS id) {

  int status;
  status=dHIPS_Clean(id);
  return(status);  
};

INTS HIPS_Clean(hips_t< complex<double> > *parameters, INTS id) {

  int status;
  status=zHIPS_Clean(id);
  return(status);  
};

INTS HIPS_Finalize(hips_t<double> *parameters) {

  int status;
  status=dHIPS_Finalize();
  return(status);  
};

INTS HIPS_Finalize(hips_t< complex<double> > *parameters) {

  int status;
  status=zHIPS_Finalize();
  return(status);  
};

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HIPSget_IDs(hips_t<double> *parameters)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id=0;
  while (DHIPS_SOLVER[id]!=-1) id++;
  DHIPS_SOLVER[id]=id;
  
  return(id);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HIPSget_IDs(hips_t< complex<double> > *parameters)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int id=0;
  while (ZHIPS_SOLVER[id]!=-1) id++;
  ZHIPS_SOLVER[id]=id;
  
  return(id);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HIPS_Initialize(hips_t<double> *parameters, int nmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=dHIPS_Initialize(nmax);
  dHIPS_ExitOnError(status);
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HIPS_Initialize(hips_t< complex<double> > *parameters, int nmax)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n;
  int status;
  status=zHIPS_Initialize(nmax);
  zHIPS_ExitOnError(status);
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HIPSinitialise_IDs(hips_t<double> *parameters)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, n;
  
  if(DHIPS_SOLVER==0) {
    DHIPS_SOLVER=(int *)  malloc(100*sizeof(int));
    for(n=0;n<100; n++) DHIPS_SOLVER[n]=-1;
    }
    
  if (HIPS_SOLVERNB == 0) {
    int idnbrmax = 100;
    status=HIPS_Initialize(parameters, idnbrmax);
    }
    
  HIPS_SOLVERNB++;
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HIPSinitialise_IDs(hips_t< complex<double> > *parameters)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, n;
  
  if(ZHIPS_SOLVER==0) {
    ZHIPS_SOLVER=(int *) malloc(100*sizeof(int));
    for(n=0;n<100; n++) ZHIPS_SOLVER[n]=-1;
    }
    
  if (HIPS_SOLVERNBZ == 0) {
    int idnbrmax = 100;
    status=HIPS_Initialize(parameters, idnbrmax);
    }
    
  HIPS_SOLVERNBZ++;
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HIPSfree_IDs(hips_t<double> *parameters, int id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  DHIPS_SOLVER[id]=-1;
  HIPS_SOLVERNB--;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 Free HIPS internal structure

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(HIPS_SOLVERNB==0) {
    status = HIPS_Finalize(parameters);
    HIPS_ExitOnError(parameters, status);
    }
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int HIPSfree_IDs(hips_t< complex<double> > *parameters, int id)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  ZHIPS_SOLVER[id]=-1;
  HIPS_SOLVERNBZ--;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 Free HIPS internal structure

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(HIPS_SOLVERNBZ==0) {
    status = HIPS_Finalize(parameters);
    HIPS_ExitOnError(parameters, status);
    }
  
  return(0);
}

#endif

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  hips interface
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int hips_factorize_template(SOLVER_t<T> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status,ierr;
#ifdef HIPS
#ifdef HAVE_MPI
  hips_t<T> *parameters;
  INTS  sym_pattern, sym_matrix;
  INTS id_hips, idnbr, i, j;
  INTS *unknownlist, *hipsnodelist;
  T *x, *rhsloc;
  INTS proc_id, n, ln;
  INTL *ia, nnz;
  INTS *ja;
  T *a;
  INTS domsize, nprocs;
  INTS pbegin, pend;
  INTS *mapptr,  *mapptr2,*mapp,*mapp2;
  INTS p;
  triplet_t<T> *Mat;
  bool talk;
  int root=0;

  MPI_Comm communicator;
  int nglob, nnzglob, tmp, bcl, idi, info;
  double tol;
  
  parameters = (hips_t<T> *) solver->parameters; 
  
/*------------------------------------------------------------------------------
  handle reduced cpu set case (small system do) */
  if(solver->rank==-1) return(0);
    
  id_hips=parameters->id_hips;
  communicator=solver->communicator;
  
  MPI_Comm_rank(communicator, &proc_id);
  MPI_Comm_size(communicator, &nprocs);
  
  status = MPI_Barrier(communicator);
  
  talk=(debug or (verbose==1 && proc_id==0));
  
  Mat=solver->COOtriplet;
  
  /* nglob */
  tmp=Mat->nrows;
    
  MPI_Allreduce( &tmp, &nglob, 1, MPI_INT, MPI_SUM, communicator );
  
  if (talk) printf("%s cpu=%d (over %d) : id_hips=%d nloc=%d nglob=%d nnz=%d\n", __FUNCTION__, proc_id, nprocs, id_hips, tmp, nglob, Mat->nnz);
  
/*------------------------------------------------------------------------------
  parameter domsize is an argument of testHIPS.exe */
  domsize = nglob/nprocs/2;

  if (talk) printf("%s cpu=%d : domsize=%d\n",__FUNCTION__,proc_id,domsize);

  HIPS_SetOptionINT(parameters, id_hips, HIPS_PARTITION_TYPE, 0); 
  HIPS_SetOptionINT(parameters, id_hips, HIPS_DOMSIZE, domsize);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  OPTIONS
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   ierr= HIPS_SetDefaultOptions (id_hips, HIPS_HYBRID); 
  ierr=HIPS_SetDefaultOptions(parameters, id_hips, HIPS_ITERATIVE);
  
  tol=1e-12;
  HIPS_SetOptionREAL(parameters, id_hips,HIPS_DROPTOL0,tol);
  HIPS_SetOptionREAL(parameters, id_hips,HIPS_DROPTOL1,tol);
  
//   tol=1e-15;
//   HIPS_SetOptionREAL(parameters, id_hips,HIPS_PREC,tol);  // HIPS TESTING

  if(verbose==1)
    info=4;
  else
    info=0;
  HIPS_SetOptionINT(parameters, id_hips, HIPS_VERBOSE, info);
  
  sym_matrix=0;
  HIPS_SetOptionINT(parameters, id_hips, HIPS_SYMMETRIC, sym_matrix);
/*------------------------------------------------------------------------------
  C : numbering starts from 0 but here starts from 1 */
  HIPS_SetOptionINT(parameters, id_hips, HIPS_FORTRAN_NUMBERING, 0);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  ENTER THE GRAPH
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  ierr = HIPS_GraphBegin(parameters, id_hips, nglob, Mat->nnz);
  HIPS_ExitOnError(parameters, ierr);
  
  for (bcl=0;bcl<Mat->nnz;bcl++) {
    ierr=HIPS_GraphEdge(parameters, id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1);
    }
  ierr = HIPS_GraphEnd(parameters, id_hips);
  HIPS_ExitOnError(parameters, ierr);
  
  status = MPI_Barrier(communicator);
  
  if (talk) printf("%s cpu=%d : GraphEnd OK...\n",__FUNCTION__,proc_id);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
   ENTER A USER PARTITION
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  mapptr = (INTS *)malloc(sizeof(INTS)*(nprocs+1));
  mapptr2= (INTS *)malloc(sizeof(INTS)*(nprocs+1));
  
  mapp =   (INTS *)malloc(sizeof(INTS)*nglob);
  mapp2=   (INTS *)malloc(sizeof(INTS)*nglob);

  for (bcl=0;bcl<nprocs;bcl++)  mapptr2[bcl]=0;
  for (bcl=0;bcl<nglob;bcl++)   mapp2[bcl]  =0;
  for (bcl=0;bcl<nprocs;bcl++)  mapptr[bcl] =0;
  for (bcl=0;bcl<nglob;bcl++)   mapp[bcl]   =0;
  
  for (bcl=proc_id+1;bcl<nprocs+1;bcl++) mapptr2[bcl]=Mat->nrows;

  MPI_Allreduce(mapptr2, mapptr, nprocs+1, MPI_INT, MPI_SUM, communicator );
 
/*------------------------------------------------------------------------------
  handle empty partition case */
  if(Mat->nnz!=0) {
    idi=0;
    mapp2[idi+mapptr[proc_id]] =  Mat->i[0]-1;
    for (bcl=1; bcl < Mat->nnz; bcl++) {
      if ( mapp2[idi+mapptr[proc_id]] != (Mat->i[bcl]-1) ) {
        idi++;
        mapp2[idi+mapptr[proc_id]] = Mat->i[bcl]-1;
        }
      }
    }
    
/*------------------------------------------------------------------------------
  Only the master processor (0 by default) needs to enter the partition */
//   MPI_Allreduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, communicator );
  MPI_Reduce( mapp2, mapp, nglob, MPI_INT, MPI_SUM, root, communicator );
            
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Set the unknownlist 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  handle empty partition case */
  if(Mat->nnz!=0) {
    unknownlist=(INTS *) malloc(sizeof(INTS)*(Mat->nrows));
    unknownlist[0]=Mat->i[0]-1;
    idi=0;
    for (bcl=1; bcl < Mat->nnz; bcl++) {
      if ( unknownlist[idi] != (Mat->i[bcl]-1) ) {
        idi++;
        unknownlist[idi] = Mat->i[bcl]-1;
        }
      }
    }
  else unknownlist=(INTS *) malloc(sizeof(INTS)*(1));

  if (talk) printf("%s cpu %d : ndofs=%d (over %d)\n",__FUNCTION__,proc_id,idi+1,Mat->nrows);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  Enter Partition 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if (proc_id ==0) {
    ierr=HIPS_SetPartition(parameters, id_hips, nprocs, mapptr, mapp);
    HIPS_ExitOnError(parameters, ierr);
    }
  
  status = MPI_Barrier(communicator);
  if (talk) printf("%s cpu %d : HIPS_SetPartition done\n",__FUNCTION__,proc_id);
  
  free(mapptr2);
  free(mapptr);
  free(mapp2);
  free(mapp);
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  GET THE LOCAL UNKNOWN LIST
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   ierr = HIPS_GetLocalUnknownNbr(id_hips, &ln);
//   HIPS_ExitOnError(ierr);
//   hipsnodelist = (INTS *)malloc(sizeof(INTS)*ln);
//   ierr = HIPS_GetLocalUnknownList(id_hips, hipsnodelist);
//   HIPS_ExitOnError(ierr);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  ENTER THE MATRIX COEFFICIENT
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  /** The processor enter the rows pbegin to pend-1 (naive partition) **/
  /** In the HIPS_ASSEMBLY_FOOL mode any processor can enter any coefficient :
      nevertheless it is better to enter as much as possible coefficients
      that are in the local matrix of the processor **/

  ierr=HIPS_AssemblyBegin(solver, id_hips, Mat->nnz, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_OVW, dHIPS_ASSEMBLY_FOOL, sym_matrix);
  HIPS_ExitOnError(parameters, ierr);
  
  if (talk) printf("%s cpu %d : AssemblyBegin OK...\n",__FUNCTION__,proc_id);

  for(bcl=0;bcl<Mat->nnz;bcl++) {
    ierr = HIPS_AssemblySetValue(id_hips, Mat->i[bcl]-1, Mat->j[bcl]-1, Mat->x[bcl]);
    HIPS_ExitOnError(parameters, ierr);
    }

  ierr = HIPS_AssemblyEnd(solver, id_hips);
  HIPS_ExitOnError(parameters, ierr);
   
  parameters->unknownlist = unknownlist;
  parameters->n   = Mat->nrows;
  parameters->nnz = Mat->nnz;

  if (talk) printf("%s cpu %d : done...\n",__FUNCTION__,proc_id);

  status = MPI_Barrier(communicator);
  
  ierr=0;
  status=0;
  return(status);  
  
#endif 
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hips_factorize(SOLVER_t<double> *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=hips_factorize_template(solver, verbose, debug);
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hips_factorize(SOLVER_t< complex<double> > *solver, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  debug=true;  // HIPS TESTING
  verbose=1;   // HIPS TESTING
  status=hips_factorize_template(solver, verbose, debug);
  return(status); 
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int hips_solve_template(SOLVER_t<T> *solver, T *RHS, int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
#ifdef HIPS
#ifdef HAVE_MPI
  hips_t<T> *parameters;
  char *solver_name;
  int proc_id, nprocs, ierr, status, n, nnz,bcl, nglob;
  int root_cpu=0;
  double cput1, cput2, cput3;
  INTS id_hips;
  INTS *unknownlist;
  MPI_Comm communicator;
  
/*------------------------------------------------------------------------------
  handle reduced cpu set case (small system do) */
  if(solver->rank==-1) return(0);

  parameters = (hips_t<T> *) solver->parameters;
  id_hips =parameters->id_hips;
  n=parameters->n;

  communicator=solver->communicator;
  
  cput1 = MPI_Wtime();
  ierr  = MPI_Comm_rank(communicator, &proc_id);
  ierr  = MPI_Comm_size(communicator, &nprocs);
  
  if(debug) printf("%s cpu=%d (over %d) : id_hips=%d\n", __FUNCTION__, proc_id, nprocs, id_hips);

  nnz=parameters->nnz;
  unknownlist = parameters->unknownlist;
    

  if(debug) printf("%s cpu=%d (over %d) : set RHS\n", __FUNCTION__, proc_id, nprocs, id_hips);
  
  int mode=solver->RHS_distribution; 
  switch(mode) {
    case SOLVER_CENTRALIZED:
/*------------------------------------------------------------------------------
      root_cpu is the global rhs passing processor */
      ierr = HIPS_SetGlobalRHS(id_hips, RHS, root_cpu, HIPS_ASSEMBLY_OVW);
      HIPS_ExitOnError(parameters, ierr);
      break;
    case SOLVER_DISTRIBUTED:
      ierr = HIPS_SetRHS (id_hips, n, unknownlist, RHS,  HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_FOOL);
//       ierr = HIPS_SetRHS (id_hips, n, unknownlist, RHS,  HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_OVW, HIPS_ASSEMBLY_RESPECT);
      HIPS_ExitOnError(parameters, ierr);
      break;
    }

  cput1 = MPI_Wtime();
  
  if(debug) printf("%s cpu=%d (over %d) : get solution\n", __FUNCTION__, proc_id, nprocs, id_hips);
  
  mode=solver->solution_distribution; 
  switch(mode) {
    case SOLVER_CENTRALIZED:
      ierr = HIPS_GetGlobalSolution(id_hips,RHS,-1);
      HIPS_ExitOnError(parameters, ierr);
      break;
    case SOLVER_DISTRIBUTED:
      ierr = HIPS_GetSolution(id_hips,n,unknownlist, RHS, HIPS_ASSEMBLY_FOOL);
      HIPS_ExitOnError(parameters, ierr);
      break;      
    }
  
//   HIPS_SetOptionINT (id_hips,HIPS_DISABLE_PRECOND, 1);

  cput2 = MPI_Wtime();
  cput3 = cput2-cput1;
  
#ifdef VERBOSE  
  if (myid == 0) {
    printf("  HIPS Resolution time =%lf \n",cput3);
  }
#endif 

  if(debug) printf("%s cpu=%d (over %d) : done...\n", __FUNCTION__, proc_id, nprocs, id_hips);
  
  if(ierr==1) ierr=0;
  else {
    printf("%s cpu=%d (over %d) : failed...\n", __FUNCTION__, proc_id, nprocs, id_hips);
    }
  return(ierr);
  
#endif
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hips_solve(SOLVER_t<double> *solver, double *RHS, double *x,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=hips_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hips_solve(SOLVER_t< complex<double> > *solver, complex<double> *RHS, complex<double> *x,int transposed, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  debug=true;  // HIPS TESTING
  status=hips_solve_template(solver, RHS, transposed, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int hips_initialize_template(SOLVER_t<T> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm communicator;
  int rank;
    
#ifdef HIPS
  hips_t<T> *parameters=new hips_t<T>;
  solver.parameters = parameters;
  
  status=HIPSinitialise_IDs(parameters);
  parameters->unknownlist=0;
  
  parameters->id_hips=HIPSget_IDs(parameters); /** id of the linear system **/
  
  HIPS_SOLVERNB++;
  if(solver.rank!=-1) {
    status = HIPS_SetCommunicator(parameters, parameters->id_hips, solver.communicator);
    HIPS_ExitOnError(parameters, status);
    }
    
  solver.targeted_packing    = MATRIX_PACKED;
  solver.targeted_ordering   = MATRIX_ROW_MAJOR;
  solver.targeted_addressing = MATRIX_TRIPLET;
  solver.targeted_numbering  = MATRIX_F_NUMBERING;
  
  solver.initialize = hips_initialize;
  solver.terminate  = hips_terminate;
  solver.factorize  = hips_factorize;
  solver.solve      = hips_solve;
  
  solver.id=SOLVER_ID_HIPS;
  
#else
    printf("Please compile poc-solvers library with -DHIPS \n");
    exit(-1);
#endif
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hips_initialize(SOLVER_t<double> & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=hips_initialize_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hips_initialize(SOLVER_t< complex<double> > & solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=hips_initialize_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int hips_terminate_template(SOLVER_t<T> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  MPI_Comm communicator;
  
#ifdef HIPS
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

 Free HIPS internal structure for problem id

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  hips_t<T> *parameters = (hips_t<T> *) solver->parameters;
  
  int id_hips=parameters->id_hips;
  
  status = HIPS_Clean(parameters, id_hips);
  HIPS_ExitOnError(parameters, status);
  
  HIPSfree_IDs(parameters, id_hips);
  
  if(status=dHIPS_SUCCESS) status=0;
  
  parameters->destroy();
  
  return(status);
#else
  return(status);
#endif 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hips_terminate(SOLVER_t<double> *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=hips_terminate_template(solver, verbose, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hips_terminate(SOLVER_t< complex<double> > *solver, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  
  status=hips_terminate_template(solver, verbose, debug);
  
  return(status);
}



/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  solver interface compatibility with older versions
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

#if OLD_SOLVER_VERSION
// if activated, use historical (1.0) code (still available in the c-code branch)
// not maintained, kept only for testing purposes, aimed to be removed
#else
// if not, provide historical (1.0) interface to new code (2.0)
// designed for smooth transition purposes, aimed to be removed

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorize_hips(solver_t *slv, triplet *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);

  slv->Mat = Mat; 
      
  solver.d_import2packed(slv);

  status=hips_factorize(&solver, verbose, debug);
  
  solver.d_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int factorizez_hips(solver_t *slv, tripletz *Mat, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  slv->Mat = Mat;
  
  solver.z_import2packed(slv);

  status=hips_factorize(&solver, verbose, debug);
  
  solver.z_export_packed(slv);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solve_hips(solver_t *slv, double *RHS, int transposed, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int verbose=0;
  
  SOLVER_t<double> solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.d_import2packed(slv);

  status=hips_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solvez_hips(solver_t *slv, complex<double> *RHS, int transposed, int verbose, bool debug)
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  SOLVER_t< complex<double> > solver;
 
  if(verbose==1) printf("enter %s\n",__FUNCTION__);
      
  solver.z_import2packed(slv);

  status=hips_solve(&solver, RHS, 0, transposed, verbose, debug);
  
  return(status); 
}

#endif

