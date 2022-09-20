
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

//extern "C" { 
#include "solverlib.h"
//}


//#include "revision.def"
#include "init_config.hpp"

#ifndef TUGO_COMMONS_H
#define TUGO_COMMONS_H

#define TUGO_VERSION (TUGO_MAJOR_VERSION*10000+TUGO_MINOR_VERSION*100+TUGO_REVISION_VERSION)

#ifdef MAIN_SOURCE
#   define PREFIX
#   define INITIALVALUE(x) =x
#else
#   define PREFIX extern 
#   define INITIALVALUE(x)
#endif

/**-----------------------------------------------------------------------------
  max number of subcycling levels */
#define NLVL 4
/**-----------------------------------------------------------------------------
  max number of dynamically stored time frames */
#define NTFMAX 3

#define COLD_START 0
#define HOT_START  1
#define EXTERNAL   2

#define CNES_TIME  0
#define MODEL_TIME 1

#define KH_CONSTANT        0
#define KH_PROPORTIONAL    1
#define KH_SMAGORINSKY     2
#define KH_SMAGORINSKY_2D  3
#define KH_UPWINDING       4
#define KH_NONE            5

#define KV_CONSTANT        0
#define KV_YAMADA          1
#define KV_GALPERIN        2
#define KV_GASPAR          3
#define KV_TKE             4
#define KV_KOMEGA          5
#define KV_HOMOGENEOUS     6 

#define KV_STRG_CONSTANT    "CONSTANT"
#define KV_STRG_YAMADA      "YAMADA"
#define KV_STRG_GALPERIN    "GALPERIN"
#define KV_STRG_GASPAR      "GASPAR"
#define KV_STRG_TKE         "TKE"
#define KV_STRG_KOMEGA      "KOMEGA"
#define KV_STRG_HOMOGENEOUS "HOMOGENEOUS"

#define GWE_CLX            0
#define GWE_XXX            1
#define WE_EXPLICIT        2
#define WE_IMPLICIT        3
#define FV_IMPLICIT        4

#define FINITE_ELEMENTS    0
#define FINITE_VOLUMES     1

#define FORWARD            0
#define LEAPFROG           1

#define BOUSSINESQ         0
#define COMPRESSIBLE       1

#define NODE_WISE          0
#define ELEMENT_WISE       1

#define PRIMITIVE          0
#define CONSERVATIVE       1
#define LAGRANGIAN         2

#define UNDEFINED_MODE     0
#define BAROTROPIC_MODE    1
#define BAROCLINIC_MODE    2
#define SW2D1D_MODE        3

#define WW3_SG_NETCDF 0
#define WW3_UG_NETCDF 1
#define WW3_SG_NATIVE 2
#define WW3_UG_NATIVE 3

#define SEQUENTIAL_COMPUTING  0
#define PARALLEL_COMPUTING    1

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : MANDATORY !!!

  Note: MPI PARALLEL COMPUTING

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

#define EXTERNAL_Z_BOUNDARY -1
#define INTERNAL_Z_BOUNDARY -2

#define EXTERNAL_U_BOUNDARY -1
#define INTERNAL_U_BOUNDARY -2

#define OBC_CLAMPED        1
#define OBC_FLATHER        2
#define OBC_MIXTE          7

#define SEQUENTIAL_OUTPUT  0
#define GLOBALIZED_OUTPUT  1
#define PARTITIONED_OUTPUT 2

#define BOUNDARY_LOCAL       0
#define BOUNDARY_VARIATIONAL 1

#define SEPARATOR_1 "\nxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx\n"
#define SEPARATOR_2 "\n-------------------------------------------------------------------------------\n"

/**--------------------------------------------------------------------------
Global variables */
PREFIX int TUGOm_initialized;
PREFIX int gIteration;

PREFIX int we_type, me_type, gFormulation, gMomentumIntegration;

PREFIX int u2D_discretisation,z2D_discretisation,g2D_discretisation;
PREFIX int u3D_discretisation,z3D_discretisation,w3D_discretisation,t3D_discretisation;

PREFIX int gDiscretisation2D,gDiscretisation3D;
PREFIX int gSequentialPaire2D;

PREFIX int advection2D_mode; 
PREFIX int AD_discretisation;

PREFIX int sproduct2D;
PREFIX int turbulence_type;
PREFIX int compressibility;

PREFIX int AutomaticInterval;

PREFIX tugo_cfg_C *tugo_cfg;
PREFIX char *parameter_file;

PREFIX int gNTF;

PREFIX int gGravityMode, gAnelasticBodyTide;

PREFIX int    gCoriolisMode;
PREFIX double gCoriolisLatitude;

PREFIX int gUseReformedAction;

/**--------------------------------------------------------------------------
gstate3D is an intermediate state vector, need to be clarified*/
PREFIX state_t      gNCP1state3D,baroclinic_delta3D;
/**--------------------------------------------------------------------------
3D state vector*/
PREFIX state_t          gstate3D[3];
PREFIX state_t          gtmp2D[3];
PREFIX action_t     gaction3D;
/**--------------------------------------------------------------------------
2D state vector, gmean2D contains time-averaged 2D state vector*/
PREFIX state2d_t      gstate2D[NLVL][3];
PREFIX parameter_t    gdata2D [NLVL];
PREFIX state2d_t      gmean2D[NLVL];
PREFIX action_t       gAction2D[NLVL][NTFMAX];

PREFIX state2d_t      gP1state2D;
PREFIX parameter_t    gP1data2D;
PREFIX int            gP1discretisation;

PREFIX climatology2D_t  climatology2D;
PREFIX state_t          climatology3D;

PREFIX double ***vertex_pressure, ***edge_pressure;

PREFIX mesh_t gFEmesh[NLVL];
PREFIX mesh_t *gMainFEmesh;

/**----------------------------------------------------------------------------
 parallel computing variables (MPI) */
// PREFIX int gParallel;
// PREFIX int gCPU_ID, gCPU_MASTER;
// PREFIX int gnCPUs;

PREFIX int gArchiveExtent;
PREFIX int gRestartExtent;

/**----------------------------------------------------------------------------
 parallel computing variables (MP) */
PREFIX int gOPENMP;
PREFIX int gOPENMP_nCPUs, gOPENMP_nCPUsMax;

/**----------------------------------------------------------------------------
 boundary specifications */
PREFIX char BelFile[1024],PeriodicFile[1024];

PREFIX z_spec  *gZspec[NLVL];
//PREFIX l_bndry *gFspec[NLVL];

//PREFIX int gNobc_F[NLVL], gNobc_Z[NLVL];
//PREFIX int gNobc_Z[NLVL];
PREFIX int NTidalData, NTidalWave, **bdy_nde;
PREFIX int *F1nd, *L2nd, *CalcCode, **NdsCds;
PREFIX int *boundary_codes;

PREFIX double *boundary_matrix[NLVL];
PREFIX int    *boundary_pivot[NLVL], boundary_hbw[NLVL];

PREFIX specification_obsolete_t gOpenBCs[NLVL], gFluxBCs[NLVL];
PREFIX specification_obsolete_t gSpectralOpenBCs, gSpectralRigidBCs;

PREFIX double Z2_ramp, Z2_time;

PREFIX sample_set_t  gSampleSet;
PREFIX section_set_t gSectionSet;
PREFIX tidaldata_t  *TidalData;

// PREFIX riverspec_t *grivers[NLVL];
// PREFIX int         gnrivers[NLVL];
PREFIX vector <river_t> gRiversInput[NLVL];
PREFIX vector <river_t> gConstraints[NLVL];

/**----------------------------------------------------------------------------
 atmospheric forcing */
/**----------------------------------------------------------------------------
 removing ancient structures */
PREFIX meteo_t  meteo1, meteo2, meteo3, meteo_mean;
PREFIX int      meteo_nvar, meteo_step;
PREFIX short   *meteo_landsea;

PREFIX int  meteo_format;
PREFIX char meteo_formatname[1024], meteo_convention[1024];
PREFIX char meteo_maskfile[1024];

/**----------------------------------------------------------------------------
 ocean waves */
PREFIX wave_t *waves1, *waves2, *waves3;
PREFIX int     waves_nvar;
PREFIX short  *waves_landsea;
PREFIX char    gWavesDirectory[1024];

PREFIX int  g_waves_format;
PREFIX char g_waves_formatname[1024], gWavesConvention[1024];
PREFIX char g_waves_maskfile[1024];

#   define HOURLY  0
#   define DAILY   1
#   define MONTHLY 2
#   define WEEKLY  4
#   define YEARLY  3

PREFIX int meteo_type;

/*--------------------------------------------------------------------------
       element based arrays
--------------------------------------------------------------------------*/

/*--------------------------------------------------------------------------
assembled matrix A (L.H.S.) and (R.H.S.) vector B of Ax = B
--------------------------------------------------------------------------*/

PREFIX double *Zb;

/* UMFPACK Solver */
// PREFIX void *A_symbolic[NLVL];
// PREFIX void *A_numeric[NLVL];
// PREFIX ordering_t  A_ordering[NLVL];

#   ifdef UMFPACK
PREFIX double A_Info[UMFPACK_INFO][NLVL], A_Control[UMFPACK_CONTROL][NLVL];
#   endif
       // UMFPACK

#   ifdef PETSEC
/* Iterativ Solver */
PREFIX Mat Aiter[NLVL];
PREFIX KSP ksp[NLVL];
PREFIX Vec xiter[NLVL];

#   endif
       // PETSEC

       //solvlib 
//PREFIX triplet WE_triplet[NLVL];
//PREFIX solver *WE_solver[NLVL];
PREFIX double **invAu[NLVL], **invAv[NLVL];

//#   ifdef WEDBLE

PREFIX double *B;
PREFIX double *X;
PREFIX double **smb1, **smb2;
PREFIX double *gRHS_1, *gRHS_2;
PREFIX double *AGB;
PREFIX double *subAGB, *subB;

PREFIX double *A_banded[NLVL];
//PREFIX double *A_packed[NLVL];

PREFIX hypermatrix_t gWEMatrix[NLVL];

// #   else
// 
// PREFIX float *B;
// PREFIX float **smb1, **smb2;
// PREFIX float *AGB;
// PREFIX float *subAGB, *subB;
// 
// PREFIX float *A_banded[NLVL];
// //PREFIX float *A_packed[NLVL];
// #   endif

PREFIX int *pivot[NLVL];

PREFIX matrix2x2_t  *gDragMatrix[NLVL];
PREFIX zmatrix2x2_t *gFricMatrix;

/*--------------------------------------------------------------------------
buffer vectors needed for model gradient, divergence etc... computation
--------------------------------------------------------------------------*/
PREFIX double *gZx, *gZy, *divU, *Hcor;
PREFIX int *gZn;

/**--------------------------------------------------------------------------
harmonic analysis*/
PREFIX double gHarmonicAcquisitionInterval, gHarmonicAcquisitionStart, gHarmonicFirstAcquisitionStart;
PREFIX double gHarmonicAnalysisInterval;
PREFIX int    gHarmonic_StructuredArchive, gHarmonic_LGP1Archive, gHarmonic_NETCDFArchive;
PREFIX int    gHarmonic_RecycleSequential;

PREFIX double *gHarmonicMatrix/* STS:, **bharm[3] already declared */;
PREFIX hrhs_t  gHarmonicRHS;
PREFIX int nharm;

/**--------------------------------------------------------------------------
Subcycling*/
PREFIX float *gCFL;
PREFIX int *history[NLVL+1], *gMaxlevel, gMaxdepth;
PREFIX int *gLocked, gLockLimit;
PREFIX double gCFLratio;
PREFIX int gNsuspects;
PREFIX int nsubcycles;
PREFIX int smooth_instabilities;
PREFIX int persistence;

PREFIX int    modelStart;
PREFIX double gMinDepth,gRelativeMinDepth,gDryingMinDepth;;

/*--------------------------------------------------------------------------
file names and other names*/
PREFIX char GridName[1024], GridNameRoot[1024], gMeshFileFormat[10];
PREFIX Cgrid_t gCgrid;

PREFIX int  gMeshFileType;

PREFIX char runid[1024], run_list_file[1024];

PREFIX char UGarchive_rootname[1024], UGforcing_rootname[1024], *gOutputPath, g;

PREFIX char ContinueFile[1024];
PREFIX char HotStartFile[1024];
PREFIX char PlotInputFile[1024];
PREFIX char SectionInputFile[1024];

PREFIX char gRecordFile[1024], trackingFile[1024];

PREFIX char ForcingD[1024];
PREFIX char LSA_directory[1024], LSA_convention[1024];
PREFIX char SAP_directory[1024], SAP_convention[1024];

PREFIX char TideDataFile[1024];

PREFIX char TopoFile[1024], SlopeFile[1024], ZminFile[1024], RugosityFile[1024];
PREFIX char BVFrequencyFile[1024], FrictionFile[1024], BackgroundFile[1024];
PREFIX float BVFrequencyDefault;
PREFIX char WaveFile[1024], IceCoverFile[1024],IceElasticityFile[1024];

PREFIX string TEnergyFile, TrackingFile;

PREFIX char tugo_version[256];

/*--------------------------------------------------------------------------
polygones*/
PREFIX char *subdomain;
PREFIX plg_t *polygones;
PREFIX int npolygones;

/*--------------------------------------------------------------------------
Time variables:

* dT is the run time step
* t is the time in seconds elapsed from t_reference
* t_end  time stop for the run in seconds elapsed from t_reference,
  computed from the start time + run duration

--------------------------------------------------------------------------*/

PREFIX double t_start, t_end, t_current; 
PREFIX date_t t_reference;  /// the model start date and time, read from model config
PREFIX double t_nodal;
//PREFIX double dT, subdT;
PREFIX double deltaT[NLVL],baroclinic_dT;

/*--------------------------------------------------------------------------
Tidals wave lists : spectrum is a global list, from Darwin/Schureman
                    WaveList is the actual list used in the run for 
                    forcing and boundary conditions
--------------------------------------------------------------------------*/

PREFIX spectrum_t WaveList, AdmittanceList, AnalysisList, EquilibriumList;
PREFIX spectrum_t reference_spectrum;
PREFIX char atlas_directory[1024], atlas_convention[1024];

/*--------------------------------------------------------------------------
t_start: time origin for the run in seconds elapsed from t_reference
           --> when cold restart, it is 0
           --> when hot restart,  it is taken from the continuation file
--------------------------------------------------------------------------*/

PREFIX double RunDuration;
/*PREFIX double t_start;*/
PREFIX double ContFInt;
PREFIX int w_need_update;

PREFIX double PlotStart, PlotInterval;

/**--------------------------------------------------------------------------
  Archives
---------------------------------------------------------------------------*/
PREFIX double archive_start, archive_sampling; 
PREFIX double archive_scale1[3], archive_scale2[3];
PREFIX char   archive_formatname[1024], archive_convention[1024];
PREFIX int    archive_format;
PREFIX int    archive_type, archive_count, archive_update;
PREFIX int    gArchivesLevel INITIALVALUE(0), gArchivesMaxFrame INITIALVALUE(-1);

/**--------------------------------------------------------------------------
  Restarts
---------------------------------------------------------------------------*/
PREFIX char hotstart_formatname[1024], hotstart_convention[1024];
PREFIX int hotstart_format;

PREFIX char restart_formatname[1024], restart_convention[1024];
PREFIX int restart_format,restart_format3D;

PREFIX char solver_name[1024];
//PREFIX int solver_type;

PREFIX int gLinearSolverID, gSubLinearSolverID, gLinearSolverMode, gSpectralSolver;

PREFIX int plot_update;
PREFIX double EnergyStart, EnergyInterval;
PREFIX double compatibility_scale;

/*--------------------------------------------------------------------------
Physical constants
--------------------------------------------------------------------------*/
/* background velocity = U0 instead of V0, which is the greenwich argument */
PREFIX double gCd_default, gCd3D_default, gCd_minimum, gCd3D_minimum, gz0_default, gz03D_default, gR_default, gR3D_default;
PREFIX int    gBottomFrictionType;
PREFIX double Clinear,Cquadratic,Ckarman;
PREFIX double gCm, tau0, gTheta, gU0, minKh, minKl, minKv, momentum2D_UpwindingCoefficient;
PREFIX double tracers2D_minKh,tracers2D_UpwindingCoefficient;
PREFIX double tracers3D_minKh,tracers3D_minKv;
PREFIX int    Kv_mode,Kh_mode,Kh2D_mode,Kh2Dtracer_mode;
PREFIX double wdc1, wdc2;
PREFIX double gwdHmin, gwdHmax;
PREFIX int    wdAlgorithm;
PREFIX double Tmin,Tmax,Smin,Smax;
PREFIX double MaxdV, MaxdZ, MaxdT, MindZdT;

/*--------------------------------------------------------------------------
Numerical integration 
--------------------------------------------------------------------------*/
PREFIX double gauss_x[16], gauss_y[16], gauss_w[16];
PREFIX int gauss_n;

/*--------------------------------------------------------------------------
constant factors and utility arrays
--------------------------------------------------------------------------*/

PREFIX double **Au, **Av, *Bu, *Bv;
//PREFIX double **invAu, **invAv;
PREFIX double *Bpartial_u, *Bpartial_v;

PREFIX double **grad_opx, **grad_opy;
PREFIX double *Pax, *Pay, *Pex, *Pey, *tlx, *tly, *buf8;
PREFIX float **plt_arr, SpinUp, AttenuationFactor;

PREFIX int surface_pressure;
PREFIX int tidal_forcing, astronomic_forcing, LSA_forcing, SpectralAtmosphericPressure;
PREFIX int wind_forcing, pressure_forcing;
PREFIX int ocean_loading;
PREFIX int wave_drag, shear_drag, ice_drag;
PREFIX int gFullTideAvailable, gBoundary_FullTideAvailable, nodal_corrections;
PREFIX int mean_wind;
PREFIX int bottom_forcing;
PREFIX int g_oceanwaves_forcing;

PREFIX astro_angles_t gAstronomicAngles;/* NEEDED BY STS */

PREFIX double gAsselinC1,gAsselinC2;

PREFIX int tidal_bc, ib_bc, spongious_bc;
PREFIX int ElevationOBC_mode, ElevationOBC_type, ElevationOBC_ActualType;

PREFIX int need_tides, need_meteo;
PREFIX int admittance, equilibrium, harmonic, energetics, gSpectralEnergy, checks, sampling,
           runtime_save,compute_stream,dynamic_mode; 
PREFIX int variable_icecover,variable_iceelasticity;
PREFIX int need_gradient;
PREFIX int mean_3Dpressure;

PREFIX int drifters, tracers;
PREFIX char drifters_file[1024],tracers_file [1024];

PREFIX int external_obc, external_state, external_transport;
PREFIX char external_path[1024], external_root[1024];

PREFIX char quaker_convention[1024], quaker_directory[1024];

PREFIX int    regular_output;
PREFIX char   gStructuredZone[1024];
PREFIX double regular_start, regular_sampling; 

PREFIX int LGP1_updated;

#   define UNDEFINED -1
#   define LINEAR     0
#   define QUADRATIC  1
#   define CUBIC      2
#   define SPLINE     4


#   define LINEAR_FRICTION    0
#   define QUADRATIC_FRICTION 1
#   define KARMAN_FRICTION    2

#   define VELOCITY    0
#   define TRANSPORT   1

#define RESTART_FORMAT_ASCII  0
#define RESTART_FORMAT_NETCDF 1
#define RESTART_FORMAT_BINARY 2

PREFIX int friction_mode,momentum2d_mode;
PREFIX int time_interpolation, obc_time_interpolation, g_oceanwaves_time_interpolation;

PREFIX FILE *in4;       /*  input file for control parameters        */
PREFIX FILE *in5;       /*  input file for continuation              */
PREFIX FILE *in6;       /*  .bel file (boundary description)         */
PREFIX FILE *echo;      /*   Record File                             */
PREFIX FILE *out4;      /*   output file for total energy (RH5)      */
PREFIX FILE *out6;      /*   Point plot file                         */
//   output file for enrgy budgets and physic checks 
PREFIX FILE *out_checks, *out_budget;
PREFIX FILE *tracking;  /*  output file file for stability  tracking */

PREFIX int FlushTracking,PrintTracking;  /* flush  tracking */

PREFIX int SmoothZField,SmoothUField,SmoothAdvection,SmoothDiffusion;

PREFIX int Spectral2DRun, Spectral3DRun;

PREFIX complex<double> *gPMLattenuation_x, *gPMLattenuation_y;
PREFIX double *gBufferAttenuation;

/* Radius of the earth in m  , gravity, density, etc...  */
#include "constants.h"

#endif /* TUGO_COMMONS_H */
