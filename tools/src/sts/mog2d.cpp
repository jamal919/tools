
/**************************************************************************

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

***************************************************************************/

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "config.h"

#include "tugo-prototypes.h"
/** Fred: dry addition/change **/
#include "dry.h"
/** Fred: dry addition/change **/
#include "tmatrix.h"
#include "parallel.h"
#include "tides.h"

extern dummy_t obc_memory;
/* STS: NO sequential so NO ContinuityMatrix */

extern  int init_tsunami(state2d_t * P1_state, parameter_t  P1_data, int nndes);
extern const char *fct_runtime(char *rootname);

#ifdef PERTURB
#   include "perturb.h"
#endif
/*
#define CHKSYS
#undef CHKSYS
*/
exchange_t  gExchange;
connexion_t geConnexions;
partition_t gPartition;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void allocate_miscellaneous(mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Note:

  07/10/2009

    this section is a mine field in z discretisation si not LGP1
    Still in fixing state

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
  int n,status;
  int nndes, nelements;
  int udim;
  int zdim;

  status=paire_dimension(mesh,&zdim,&udim);

  nndes = mesh.nvtxs;

  if(need_meteo) {
  switch (time_interpolation) {
    case (LINEAR):
      status=meteo1.allocate(zdim,udim);
      if(status!=0)
        check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
      status=meteo2.allocate(zdim,udim);
      if(status!=0)
        check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
      break;

    case (QUADRATIC):
      status=meteo1.allocate(zdim,udim);
      if(status!=0)
        check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
      status=meteo2.allocate(zdim,udim);
      if(status!=0)
        check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
      status=meteo3.allocate(zdim,udim);
      if(status!=0)
        check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
      break;
    }

    status=meteo_mean.allocate(zdim,udim);
    if(status!=0)
      check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
    }

  if(g_oceanwaves_forcing){
    switch (g_oceanwaves_time_interpolation) {
      case (LINEAR):
        waves1= new wave_t;
        status=(*waves1).allocate(nndes); /// HERE !!! 
        if(status!=0)
            check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
        waves2= new wave_t;
        status=(*waves2).allocate(nndes);
        if(status!=0)
            check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
        break;
      case (QUADRATIC):
        waves1= new wave_t;
        status=(*waves1).allocate(nndes);
        if(status!=0)
            check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
        waves2= new wave_t;
        status=(*waves2).allocate(nndes);
        if(status!=0)
            check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
        waves3= new wave_t;
        status=(*waves3).allocate(nndes);
        if(status!=0)
            check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
        break;
      }
    }

#ifdef DRY
  gzNodeWet = new int[zdim];
  if (!gzNodeWet)    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  gTriangleWet = new int[gFEmesh[0].ntriangles];
  if (!gTriangleWet) check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
#endif

  B = t1Darray(B, (double) 1e+10, zdim);
  if(!B) check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  
  switch(gFEmesh[0].ntriangles) {
    case 0:
      nelements=gFEmesh[0].nquadrangles;
      break;
    default:
      nelements=gFEmesh[0].ntriangles;
      break;
    }
  smb1 = tmatrix(smb1, 1.e+10, nelements, 3);
  if(!smb1) check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  smb2 = tmatrix(smb1, 1.e+10, nelements, 3);
  if(!smb2) check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  
  gRHS_1=new double[nelements*3];
  gRHS_2=new double[nelements*3];

  Zb = new double[zdim];
  if (!Zb)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);

  grad_opx = tmatrix(grad_opx, 1.e+10, nndes, mesh.nnghm + 1);
  grad_opy = tmatrix(grad_opy, 1.e+10, nndes, mesh.nnghm + 1);

  gZx = new double[nndes];
  if (!gZx)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  gZy = new double[nndes];
  if (!gZy)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  Hcor = new double[nndes];
  if (!Hcor)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  gZn = new int[nndes];
  if (!gZn)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);

  tlx = new double[nndes];
  if (!tlx)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  tly = new double[nndes];
  if (!tly)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);

  Pax = new double[nndes];
  if (!Pax)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  Pay = new double[nndes];
  if (!Pay)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);

  Pex = new double[nndes];
  if (!Pex)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);
  Pey = new double[nndes];
  if (!Pey)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);

  divU = new double[nndes];
  if (!divU)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);

  buf8 = new double[nndes];
  if (!buf8)
    check_error(-1, "memory allocation error", __LINE__, __FILE__, 1);

  for(n = 0; n < zdim; n++) {
    Zb[n] = 0.;
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int hugo_initialize(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  date_t start, end;
  int i,k,n, status;
  int level;
  int zdim;
  int depth,uorder,zorder;
  char *RunControlFile=NULL,*RunControlDelta=NULL;
  int input_format=1;
  char msg[1024], RunControlFileCopy[1024];

  if(argc < 2)
    check_error(-1, "***** missing input file name, required as argument *****", __LINE__, __FILE__, 1);

  strcpy(tugo_version, "TUGO release 3.4 BETA, F. Lyard, LEGOS 01/09/2010");
  
  gOutputPath=NULL;

  input_format=1;
  for(n=1;n<argc;n++) {
    if(strcmp(argv[n],"-old")==0) {
      /* STS: NO old input format */
      continue;
      }
    if(strcmp(argv[n],"-delta")==0) {
      RunControlDelta = argv[n+1];
      fprintf(stderr, "***control delta file %s ****\n", RunControlDelta);
      n++;
      continue;
      }
    if(strcmp(argv[n],"--version")==0) {
      fprintf(stderr, "***version %s ****\n", tugo_version);
      continue;
      }
    if(strcmp(argv[n],"--output-path")==0) {
      /* STS: option necessary, see spectral_tugo_fortran_template_() */
      gOutputPath=strdup(argv[n+1]);
      n++;
      continue;
      }
    if(RunControlFile==NULL) {
      RunControlFile = argv[n];
      fprintf(stderr, "***control input file %s ****\n", RunControlFile);
      continue;
      }
    sprintf(msg, "command line argument not understood : %s \n", argv[n]);
    check_error(-1, msg, __LINE__, __FILE__, 1);
    }
  
#ifdef PARALLEL
#ifdef HAVE_MPI
  int rank = MPI::COMM_WORLD.Get_rank();
  int size = MPI::COMM_WORLD.Get_size();
#ifdef DEBUG_PARALLEL
  cout << "Tugo_initialize I am ....." << rank << " of " << size << endl;
#endif
#endif
#endif

  strcpy(tugo_version, "TUGO release 4.0 BETA, F. Lyard, LEGOS 19/02/2013");

/*-----------------------------------------------------------------------------
  initialize cfl*/
  gCFL=0;
  gLocked=0;

/*-----------------------------------------------------------------------------
  initialize type of element */
  me_type = -1;
  we_type = -1;
  u2D_discretisation=-1;
  z2D_discretisation=-1;

/*-----------------------------------------------------------------------------
  w computation */
  w_need_update = 0;

  obc_time_interpolation = LINEAR;

  tracers = 0;
  drifters = 0;

/*-----------------------------------------------------------------------------
  polygones for suspect imposed selection*/
  polygones = NULL;
  npolygones = 0;
  subdomain = NULL;

/*-----------------------------------------------------------------------------
  land sea mask */
  meteo_landsea=NULL;

  for(level=0;level<NLVL;level++) DryingMatrix[level]     = 0;
  for(level=0;level<NLVL;level++) ActualMatrix[level]     = 0;
  /* STS: NO sequential so NO ContinuityMatrix */

  t_start = 0.;

/**-----------------------------------------------------------------------------
  load model configuration file*/
  tugo_cfg=new tugo_cfg_C;
  echo=NULL;

  cout << "\n-----------------------------------------------------------------------------" << endl;
  cout << "Simulation configuration parameters: initialization" << endl << endl;
  int delta=0;
  if(RunControlDelta!=NULL) delta=1;
  status=configure_init(RunControlFile, input_format, tugo_cfg,delta);
  
  if(RunControlDelta!=NULL) {
    input_format=1;
    status=configure_init(RunControlDelta, input_format, tugo_cfg,0);
    }
  cout << "Simulation configuration parameters: Ok\n" << endl;

  switch (gLinearSolverMode) {
    case -1:
    case PARALLEL_COMPUTING:
#ifdef PARALLEL
#ifdef HAVE_MPI
      gLinearSolverMode=PARALLEL_COMPUTING;
/**----------------------------------------------------------------------------
      Multi-processor setup*/
      gParallel=1;
      gCPU_ID=rank;
      gnCPUs=size;
      gArchiveExtent=GLOBALIZED_OUTPUT;
      gRestartExtent=GLOBALIZED_OUTPUT;
      if (gnCPUs == 1) {
        gParallel=0;
        }
#else
/**----------------------------------------------------------------------------
      Multi-processor Debug*/
      gParallel=1;
      rank=0;
      size=4;
      gCPU_ID=rank;
      gnCPUs=size;
      gArchiveExtent=GLOBALIZED_OUTPUT;
      gRestartExtent=GLOBALIZED_OUTPUT;
      if (gnCPUs == 1) {
        gParallel=0;
        }
#endif
#endif
      break;
    }

  if(gLinearSolverMode==-1) gLinearSolverMode=SEQUENTIAL_COMPUTING;

/**-----------------------------------------------------------------------------
  STS: NO copy configuration file in runtime directory*/

/**-------------------------------------------------------------------------
  STS: NO open runtime track files */

/**-----------------------------------------------------------------------------
  load model mesh file*/
  cout << "\n-----------------------------------------------------------------------------" << endl;
  cout << "Unstructured mesh: initialization" << endl;
  if(gParallel==0) {
    gMainFEmesh=&(gFEmesh[0]);
    }
  else {
    gMainFEmesh= new mesh_t;
    *gMainFEmesh= gFEmesh[0];
    }

  switch(modelStart) {
    case HOT_START:
      /* STS: NO HOT_START */
      break;

    case COLD_START:
      status=read_mesh(gMeshFileFormat, GridNameRoot, *gMainFEmesh);
      break;
    }

  if(gParallel==0) {
//    gFEmesh[0].origin=gMainFEmesh;
    gFEmesh[0].origin=0;
    status = geometry(&gFEmesh[0]);
    }
  else {
    /* STS: NO sub-mesh */
    }
  gFEmesh[0].level=0;

#ifdef PERTURB
  read_input_user(argv);
#endif

  printf("\nCPU %d,#Number of elements   : %d \n",   gCPU_ID,gFEmesh[0].ntriangles);
  printf("CPU %d,#Number of nodes      : %d \n",     gCPU_ID,gFEmesh[0].nvtxs);
  printf("CPU %d,#Number of edges      : %d \n",     gCPU_ID,gFEmesh[0].nedges);
  printf("CPU %d,unpartitionned mesh   : %p  %p \n", gCPU_ID,gFEmesh[0].origin, gMainFEmesh);

#ifdef PARALLEL
#ifdef HAVE_MPI
  status = MPI_Barrier(MPI_COMM_WORLD);
#endif
#endif

/*-------------------------------------------------------------------------
  allocate action vectors/structures*/
  status=allocate_action2D(gFEmesh[0],0);

/*-------------------------------------------------------------------------
  allocate memory, to be refurbish*/
  allocate_miscellaneous(gFEmesh[0]);

  if(need_tides) {
    int verbose=(gCPU_ID==gCPU_MASTER);
    int debug=0;
    if(verbose) {
      cout << "\n-----------------------------------------------------------------------------" << endl;
      cout << "Tidal parameters: initialization" << endl;
      }
/**-------------------------------------------------------------------------
    initialise prescribed tidal wave list*/
    reference_spectrum = initialize_tide();//tools update
    if(strcmp(TideDataFile, "NONE") != 0)
      status=GetBdyObs(TideDataFile, compatibility_scale, verbose, debug);
    else
      status = load_wave_list(WaveFile);
/**-------------------------------------------------------------------------
    optionnally extend prescribed tidal wave list using admittance*/
    if(admittance) {
      /* STS: NO */
      }
    else {
      AdmittanceList.n = 0;
      }
/**-------------------------------------------------------------------------
    optionnally extend prescribed tidal wave list using equilibrium*/
    if(equilibrium) {
      /* STS: NO */
      }
    else {
      EquilibriumList.n = 0;
      }
    }

  if(need_tides || harmonic) {
    /* STS: Harmonic analysis is irrelevant. */
    }

/**-------------------------------------------------------------------------
  initialize 2D state vectors*/
  allocate_UGO2D(gFEmesh[0]);

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Note: MPI PARALLEL COMPUTING

  26/02/2010:

    Next instruction might be affected by parallelisation, as it processes
    the model inputs.

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
/**-------------------------------------------------------------------------
  load additional data files, using specific VERTICES database */

  zdim=gFEmesh[0].nvtxs;
  gP1data2D.lon  =new double[zdim];
  gP1data2D.lat  =new double[zdim];
  for(i = 0; i < zdim; i++) {
    gP1data2D.lon[i] = gFEmesh[0].vertices[i].lon*d2r;
    gP1data2D.lat[i] = gFEmesh[0].vertices[i].lat*d2r;
    }

  gP1state2D.vic  =new float[zdim];
  gP1state2D.vie  =new float[zdim];
  gP1state2D.N    =new float[zdim];

  gP1state2D.colocation = 1;
  
  gP1data2D.wdc1     =new float [zdim];
  gP1data2D.Nbar     =new float [zdim];
  gP1data2D.celerity =new float [zdim];
  
  gP1data2D.vwr   =new float [zdim];
  gP1data2D.h     =new double[zdim];
  gP1data2D.u0    =new double[zdim];
  gP1data2D.dhdx  =new double[zdim];
  gP1data2D.dhdy  =new double[zdim];

  gP1data2D.z0    =new float[zdim];
  gP1data2D.Cd    =new float[zdim];
  gP1data2D.rlinear     =new float[zdim];
  gP1data2D.FrictionRatio =new float[zdim];

  gP1data2D.colocation = 1;
  switch(gFEmesh[0].nature()) {
    case FE_TRIANGLE:
      gP1discretisation=LGP1;
      load_optionals(gFEmesh[0], gP1state2D, gP1data2D, LGP1, 0);
      break;
    case FE_QUADRANGLE:
      gP1discretisation=CQP1;
      load_optionals(gFEmesh[0], gP1state2D, gP1data2D, CQP1, 0);
      break;
    default:
      gP1discretisation=-1;
      load_optionals(gFEmesh[0], gP1state2D, gP1data2D, LGP1, 0);
      break;
    }
/*-------------------------------------------------------------------------
  initialize bathymetry*/
  status=topo_init(gFEmesh[0],gP1data2D);
  
  if(gParallel!=0) {
    /* STS: NO sub-mesh */
    }
/*-------------------------------------------------------------------------
  initialize internal wave drag*/
  initialize_IWdrag(gFEmesh[0], gP1data2D, gFEmesh[0].nvtxs);

/*-------------------------------------------------------------------------
  initialize horizontal diffusion*/
  /* STS: NO diffusion */

/*-------------------------------------------------------------------------
  initialize bottom friction*/
  initialize_bottomfriction2D(gFEmesh[0], gP1data2D);
  /* STS: NO ice shelves */
  
/*-------------------------------------------------------------------------
  initialize background velocity*/
  initialize_u0(gFEmesh[0], gP1data2D, LGP1, 0);

/**-------------------------------------------------------------------------
  initialize 2D state vectors*/
  initialise_UGO2D(gFEmesh[0]);

/**-------------------------------------------------------------------------
  initialize 3D state vectors*/
  if ((tugo_cfg->model->mode.face_value() != "2D") &&(tugo_cfg->model->mode.face_value() != "NO-DYNAMIC")) {
    /* STS: NO 3D */
    }
  else {
    gFEmesh[0].nlayers=1;
    gFEmesh[0].nlevels=2;
    }
  if(gParallel!=0) {
    gMainFEmesh->nlayers=gFEmesh[0].nlayers;
    gMainFEmesh->nlevels=gFEmesh[0].nlevels;
    }

/*-------------------------------------------------------------------------
  initialise variational boundary conditions matrix*/
  if(ElevationOBC_mode==BOUNDARY_VARIATIONAL) {
    status = dBoundaryMatrix(gFEmesh[0], 0, boundary_matrix[0], boundary_pivot[0], boundary_hbw[0]);
    }
  build_boundary_conditions(gFEmesh[0], gstate2D[0], gdata2D[0]);
  
  ElevationOBC_ActualType=ElevationOBC_type;
  
  if(mean_3Dpressure) {
    /* STS: NO 3D */
    }

/*-------------------------------------------------------------------------
  FES2012 */
  /* STS: NO continuation files */;

/*-------------------------------------------------------------------------
  initialize solution at Tm and T*/
  /* STS: NO start */

/*-------------------------------------------------------------------------
  init runtime state vector constraints */
  if(tugo_cfg->constraint->OnOffFlag.face_value() == (string) "TRUE") {
    /* STS: NO contraints */
    }

/*-------------------------------------------------------------------------
  initialise elevation dynamic matrix*/
  if (tugo_cfg->model->mode.face_value() != "NO-DYNAMIC") {
    /* STS: NO sequential */
    }

/*-------------------------------------------------------------------------
  initialise LGP1 horizontal gradient operator matrix*/
  VertexGradient_operator(gFEmesh[0]);

/*-------------------------------------------------------------------------
  in case not given in an external file, compute bathymetric slopes*/
  if(need_gradient == 1) {
    topo_gradient(gFEmesh[0], gP1data2D);
    }
    
/**-------------------------------------------------------------------------
  smooth the bathymetry gradient
                         - to stabilize the model
                         - to be more consistent with internal wave theory
-------------------------------------------------------------------------*/
//  smooth_gradient(gFEmesh[0], gdata2D[0], u2D_discretisation);

/*-------------------------------------------------------------------------
  no harmonic constants available yet*/
  gFullTideAvailable = 0;

/*-------------------------------------------------------------------------
  initialize tide with correct time (???) */
  if(need_tides) {
    reference_spectrum.destroy();
    reference_spectrum = initialize_tide();//tools update
    }

/*-------------------------------------------------------------------------
  load SLA input (at elevation nodes)*/
  if(LSA_forcing)
    {/* STS: NO LSA */}

/*-------------------------------------------------------------------------
  build wave drag matrix (P1 mode)*/
  /* STS: NO wave drag because no P1 */

/*-------------------------------------------------------------------------
  allocate and set default value for suspicion and stability flags */
  /* STS: NO instability */

  if(tugo_cfg->topography->OnOffFlag.actual_value == Keywords::KEY_TRUE) {
    /* STS: NO tsunami */
    }

/*-------------------------------------------------------------------------
  make the boundary stronly viscous to avoid instabilities */
  /* STS: NO instability */

//  ShowDimension(stderr);
  ShowDimension(echo);

/*------------------------------------------------------------------------
  runtime structured archiving*/
  /* STS: NO runtime archiving */

/*------------------------------------------------------------------------
  runtime unstructured archiving*/
  /* STS: NO runtime archiving */

/*-------------------------------------------------------------------------
  init the dynamic/thermodynamic ice */

#ifdef ICE
  if (use_ice){
    Feice= (ice_t *) malloc( gFEmesh[0].nvtxs*sizeof(ice_t) );
    status=init_ice(gFEmesh[0].nvtxs,gFEmesh[0].ntriangles);
    }
#endif

/*-------------------------------------------------------------------------
  init the river discharge*/
  /* STS: NO river discharge */

/*-------------------------------------------------------------------------
  PlotStart is relative to initial model time*/
  /* STS: NO sequential */
  return (0);
}
