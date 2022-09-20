
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2011

  Unstructured Ocean Grid initiative

Contributors:
 
  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  Clement Mayet      LEGOS, Toulouse, France
  Yves Soufflet      LEGOS, Toulouse, France
  Damien Allain      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

***************************************************************************/

#define MAIN_SOURCE

#ifdef PARALLEL
#ifdef HAVE_MPI
#include <mpi.h>
#endif
#endif

#include <iostream>
#include <stdarg.h> //variadic arrays

#include "config.h"

#include "functions.h"
#include "tugo-prototypes.h"
#include "dry.h"
#include "ice.h"
#include "fe.h"

extern int hugo_initialize(int argc, char *argv[]);

static int error_pending=0;

double time2D[NLVL];


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void initialize_OPENMP_tugo()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nprocs;

  gOPENMP_nCPUsMax=tugo_cfg->model->OPENMP_nCPUsMax.numerical_value<int>() ;

#ifdef HAVE_MPI
  if(gCPU_ID==gCPU_MASTER) {
    printf("OpenMP de-activated in MPI parallel mode...\n");
    gOPENMP_nCPUsMax=1;
    }
#endif

//#ifdef OMP_H
#ifdef _OPENMP 
  if(gOPENMP_nCPUsMax==-1) gOPENMP_nCPUsMax=omp_get_max_threads();
  nprocs=omp_get_max_threads();
  gOPENMP_nCPUs=MIN(nprocs,gOPENMP_nCPUsMax);
  if(gOPENMP_nCPUs>1){
    omp_set_dynamic(1);
    omp_set_num_threads(gOPENMP_nCPUs);
    gOPENMP=1;
    }
  else {
    gOPENMP=0;
    gOPENMP_nCPUs=1;
    }
#else
  gOPENMP=0;
  gOPENMP_nCPUs=1;
#endif
  
  printf("\nOPENMP activated=%d, OPENMP_nCPUsMax %d (available), OPENMP_nCPUs %d (actually used)\n\n",gOPENMP,gOPENMP_nCPUsMax,gOPENMP_nCPUs);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  const char *fct_runtime(char *rootname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string s,executable,echofile;
  int n,pos;
  FILE *echo=NULL;
  time_t call_time;
  struct tm *call_tm;
  const int dateL=32;
  char date[dateL];
  char *tmp;

  call_time=time(NULL);
  if(call_time==-1){
    __CHKERR_LINE_FILE__(errno,"time error");
    exit(-1);
    }
  call_tm=localtime(&call_time);
  if(call_tm==NULL){
    __CHKERR_LINE_FILE__(errno,"localtime error");
    exit(-1);
    }
  /**This outputs the date as yyyy-mm-dd HH:MM:SS TZ format, as shown below :
  \code /**/ // COMPILED CODE BELOW !!!
  strftime(date,dateL,"%F %T %Z",call_tm); /** \endcode */

  s.assign(rootname);

  s+="-";
  s+=date;
  
  pos=s.rfind(" ");
  while (pos!=string::npos) {
    s[pos]='_';
    pos=s.rfind(" ");
    }

  s+=".log";

  tmp=strdup(s.c_str());
  
  return(tmp);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int gravitational_constant_01(discretisation_t & descriptor)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /**----------------------------------------------------------------------------
//   http://en.wikipedia.org/wiki/Gravity_of_Earth
// -----------------------------------------------------------------------------*/
// {
//   int n, status;
//   printf("use latitude-dependant gravitational constant\n");
//   for(n=0;n<descriptor.nnodes;n++) {
//     const double s=descriptor.nodes[n].s;
//     const double c=descriptor.nodes[n].c;
//     descriptor.nodes[n].g=9.780327*(1+0.0053024*s*s-2*0.0000058*s*c*s*c);
//     }
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int gravitational_constant_00(discretisation_t & descriptor)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /**----------------------------------------------------------------------------
//   http://en.wikipedia.org/wiki/Standard_gravity
// -----------------------------------------------------------------------------*/
// {
//   int n, status;
//   printf("use uniform gravitational constant\n");
//   for(n=0;n<descriptor.nnodes;n++) {
//     descriptor.nodes[n].g=P_g;
//     }
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int gravitational_constant(discretisation_t & descriptor)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   if(tugo_cfg->physics->gravity_mode.face_value() == "CONSTANT") {
//     gGravityMode=0;
//     }
//   else if (tugo_cfg->physics->gravity_mode.face_value() == "LATITUDE-DEPENDANT") {
//     gGravityMode=1;
//     }
//   else {
//     check_error(-1, "unreckognized option for gravity mode", __LINE__, __FILE__, 1);
//     }
//   
//   switch(gGravityMode) {
//     case 0:
//       status=gravitational_constant_00(descriptor);
//       break;
//     case 1:
//       status=gravitational_constant_01(descriptor);
//       break;
//     }
//     
//   return(0);
// }


extern int cefmo_solver(mesh_t *mesh, parameter_t data);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  double start, finish, keep; 
  char *restartfile = NULL;
  int rank,size;
  double *in,*out;
  extern int fe_4maxima(void);
//   extern int solver1D(void);
//   extern int debug_vertical_mode(void);
  
  fct_echo(argc, argv);
  
/*-------------------------------------------------------------------------*//**
  set tugo_exit() to be called at return of main() or exit */
  /* STS: NO fail safe */
/*----------------------------------------------------------------------*/
  
  printf(
    "\n### COMODO STS is a frequency-domain tidal solver, copycat from its native implementation in the T-UGOm model.\n"
    "### It can operate indifferently on triangles and quadrangles grids (including structured grids).\n"
    "### The present version is limited to 2D tides, the 3D component will be added in a near future.\n"
    "\n### Contact: florent.lyard@legos.obs-mip.fr, damien.allain@legos.obs-mip.fr\n\n"
    );

#ifdef PARALLEL 
#ifdef HAVE_MPI
//  MPI::Init(argc, argv);
//  MPI::Init();
//   int required = MPI_THREAD_MULTIPLE;
  int required = MPI_THREAD_SINGLE;
  int  provided = -1;
  MPI_Init_thread(&argc, &argv, required, &provided);
  set_MPI_Init_thread_DONE(1);
//  ierr=MPI_Init(&argc, &argv);
  rank = MPI::COMM_WORLD.Get_rank();
  size = MPI::COMM_WORLD.Get_size();
  cout << "T-UGOm cpu " << rank << " of " << size << endl;
#ifdef HIPS
  FILE *fout;
/* Le nombre de solver que hips utilise. */
/**----------------------------------------------------------------------------
  MPI unsafe */
  fout=fopen("SOLVERNB","w");
  fprintf(fout,"0\n");
  fclose(fout);
#endif
#endif
#endif

  TUGOm_initialized=0;

  fe_integrale_init();

  fprintf(stderr, "^^^^^^^^^ %s -computation starts ^^^^^^^^^\n", argv[0]);
  
/**----------------------------------------------------------------------------
  Mono-processor setup*/
  gParallel=0;
  gCPU_ID=0;
  gCPU_MASTER=0;
  gnCPUs=1;
  gArchiveExtent=SEQUENTIAL_OUTPUT;
  gRestartExtent=SEQUENTIAL_OUTPUT;
    

/**----------------------------------------------------------------------------
  orthonormal basis*/
  if(gCPU_ID==0) {/* STS: status=fe_4maxima(); only for debug in tugo */}
  
  
//  status= debug_vertical_mode();

/**----------------------------------------------------------------------------
  Initialize simulation*/
  gLinearSolverID    =-1;
  gSubLinearSolverID =-1;

  gLinearSolverMode  =-1;

  gSpectralSolver    =-1;
  
  gNTF=NTFMAX;
  
  status = hugo_initialize(argc, argv);

  /* STS: NO sequential */

/**----------------------------------------------------------------------------
  number of dynamically stored time frame*/
  gNTF=3;

/**----------------------------------------------------------------------------
  OPEN MP optimisations*/
  initialize_OPENMP_tugo();
  
/**----------------------------------------------------------------------------
  Tidal spectral solver*/
  /* STS: Spectral2DRun compulsory */if(1){
    status= cefmo_solver(&gFEmesh[0],gdata2D[0]);
    }

/**----------------------------------------------------------------------------
  Sequential solver*/
  /* STS: NO sequential */

#ifdef PARALLEL
#ifdef HAVE_MPI
  MPI::Finalize();
#endif
#endif

  fprintf(stderr, "^^^^^^^^^ %s -computation complete ^^^^^^^^^\n", argv[0]);
  
  return status;
}
