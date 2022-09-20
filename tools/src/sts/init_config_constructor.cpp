
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

//  init_config_constructor.cpp 
//  
//  
//  Toutes ces fonctions sont utiles pour initialiser les classes de init_config.hpp
//  du modele TUGO
//  
//  ce fichier est indisossiable de init_config.hpp
//
//  Thierry LETELLIER LEGOS Toulouse 08Jun2006
//

#include "init_config.hpp"
#include "keywords.hpp"

model_C::model_C()
{
  rootname.implement("rootname","NONE","NONE","give the model rootname",0,"NONE","NONE");
  TimeStep.implement("time_step","NONE","NONE","the time step",0,"NONE","seconds");
  SubTimeStep.implement("sub_time_step","<default> $TimeStep/3.0","<default> $TimeStep/3.0","the sub time step",0,"NONE","second");
  OrigineTime.implement("time_origine","<default> 1950-01-01 00:00:00","<default> 1950-01-01 00:00:00","give the time origine",0,"NONE","NONE");
  MeshFormat.implement("mesh_format","<default> nei","<default> nei","give the mesh format",3,"nei nod nc","NONE");
  MeshType.implement("mesh_type","<default> TRIANGULAR","<default> TRIANGULAR","give the mesh type",5,"LINEAR TRIANGULAR QUADRANGULAR MIXED STRUCTURED","NONE");
  MeshMeta.implement("mesh_meta","<default> NONE","<default> NONE","give the mesh type",0,"NONE","NONE");
  MeshFile.implement("mesh_file","<default> $rootname.$MeshType","<default> $rootname.$MeshType","the mesh file",1,"FILE","File");

  BelFile.implement("bel_file","<default> NONE","<default> NONE","the bel file",1,"FILE","File");
  
  TopoFile.implement("topo_file","<default> NONE","<default> NONE","bathymetry file",1,"FILE","File");
  SlopeFile.implement("slope_file","<default> NONE","<default> NONE","bathymetry gradient's file",1,"FILE","File");
  ZminFile.implement("zmin_file","<default> NONE","<default> NONE","minimum elevation file",1,"FILE","File");
  RelativeMinDepth.implement("relative_minimum_depth","<default> 5.0","<default> 5.0","minimum depth relative to minimum elevation",0,"NONE","metres");
  MinDepth.implement("minimum_depth","<default> 10.0","<default> 10.0","overall minimum depth in model bathymetry",0,"NONE","metres");
  
  TopoFile.deprecated=true;
  SlopeFile.deprecated=true;
  ZminFile.deprecated=true;
  RelativeMinDepth.deprecated=true;
  MinDepth.deprecated=true;
  
  RunDuration.implement          ("run_duration","<default> 48.0","<default> 48.0","RunDuration",0,"NONE","hours");
  SpinUpDuration.implement       ("spin_up","<default> 48.0","<default> 48.0","the spin up",0,"NONE","hours");
  SpinUpMode.implement           ("spin_up_mode","<default> LINEAR","<default> LINEAR","the spin up mode",2,"LINEAR COSINE","NONE");
  RestartFile.implement          ("restart_file","<default> COLD-START","<default> COLD-START","the start condition Or filename",2,"COLD-START FILE","File");
  RestartFormat.implement        ("restart_format","<default> NETCDF","<default> NETCDF","Type of the file used for restart",2,"NETCDF CLX","File format");
  ContinuationFile.implement     ("continuation_file","<default> tugo_restart","<default> tugo_restart","Restart files rootname",0,"NONE","File");
  ContinuationFormat.implement   ("continuation_format","<default> NETCDF","<default> NETCDF","Restart file format",2,"NETCDF CLX","File format");
  ContinuationInterval.implement ("continuation_interval","<default> 10800.0","<default> 10800.0","Restart interval",0,"NONE","seconds units to be verificated");
  AutomaticInterval.implement    ("automatic_interval","<default> MONTHLY","<default> MONTHLY","Restart interval",4,"HOURLY DAILY WEEKLY MONTHLY","NONE");
  OutputPath.implement           ("output_path","<default> .","<default> .","default pathname for model outputs",1,"DIRECTORY","Directory");

  LinearSolver.implement         ("solver_type","<default> SpDOMESTIC","<default> SpDOMESTIC","Solver Type",10, "DOMESTIC LAPACK SUNPERF UMFPACK GMRES MUMPS MUMPS_SYM SpDOMESTIC HIPS PASTIX","NONE");
  SubLinearSolver.implement      ("sub_solver_type","<default> AUTO","<default> AUTO","Solver Type",10, "AUTO DOMESTIC LAPACK SUNPERF UMFPACK GMRES MUMPS MUMPS_SYM SpDOMESTIC HIPS","NONE");
  SolverMode.implement           ("solver_mode","<default> AUTO","<default> AUTO","Solver mode",3, "AUTO SEQUENTIAL PARALLEL","NONE");

  mode.implement                 ("mode","<default> 2D","<default> 2D","simulation mode",5,"NO-DYNAMIC 2D 3D-BAROTROPIC 3D-DIAGNOSTIC 3D-PROGNOSTIC","NONE");
  compressibility.implement      ("compressibility","<default> BOUSSINESQ","<default> BOUSSINESQ","compressibility mode",2,"BOUSSINESQ COMPRESSIBLE","NONE");

  WaveEqForm.implement           ("wave_equation_form","<default> P1-NCP1","<default> P1-NCP1","type of resolution",5,"P1-NCP1 P1-P1(lumped) P1-P1(exact) WE-explicit WE-implicit","NONE"); 
  WaveEqForm.deprecated=true;

  TwoDMometumEqForm.implement    ("2D_momentum_equation_form","<default> P1-NCP1","<default> P1-NCP1","type of resolution",2,"P1-NCP1 P1-P1","NONE");
  TwoDMometumEqForm.deprecated=true;

  tau0.implement("tau_0","<default> 2e-4","<default> 2e-4","weight of CE in GWE",0,"NONE","s-1");
  tau0.deprecated=true;

  theta.implement("theta","<default> 0.75","<default> 0.75","implicit/explicit coefficient",0,"NONE","units to be verificated");
  theta.deprecated=true;

  OPENMP_nCPUsMax.implement("OPENMP_nCPUsMax","<default> -1","<default> -1","Max number of CPUs allowed in OPENMP optimisation",0,"NONE","NONE");
  MetisPath.implement("MetisPath","<default> /home/softs/metis-4.0","<default> /home/softs/metis-4.0","path to METIS executables",1,"DIRECTORY","Directory");
  
  ThreeDPressure.implement("3d_pressure","<default> FALSE","<default> FALSE","FLAG",2,"FALSE TRUE","NONE");
//  ThreeDPressure.deprecated=true;

  EchoFile.implement("echo_file","<default> tugo.echo","<default> tugo.echo","the echo filename",1,"FILE","File");
  CompatibilityScale.implement("CompatibilityScale","<default> 1.0","<default> 1.0","FLAG",0,"NONE","NONE");
  CompatibilityScale.deprecated=true;

  rootname.implement_more_comments("The root correspond to BLA BLA BLA BLA BLA BLA \nBLA BLA BLA BLA \nBLA BLA BLA BLA \nBLA BLA BLA BLA \nBLA BLA BLA BLA \n" );
  MeshType.implement_more_comments("The Mesh type  correspond to BLA BLA BLA BLA BLA BLA \nBLA BLA BLA BLA \nBLA BLA BLA BLA \nBLA BLA BLA BLA \nBLA BLA BLA BLA \n" );

  OrigineTime.implement_more_comments("Time origin of model simulation \n Format YYYY-MM-DD HH:MM:SS \n" );

  RunDuration.implement_more_comments(" Simulation duration");
  SpinUpDuration.implement_more_comments(" Ramp duration for cold start");
  SpinUpMode.implement_more_comments(" shape of the initial ramp (when doing COLD-STAR)");
  MinDepth.implement_more_comments(" Overall minimum depth");
  RestartFile.implement_more_comments("input file that contains model restart data when doing hots start");
  RestartFormat.implement_more_comments(" CLX is original 2D restart format, NETCDF extends to 3D simulation");
  ContinuationFile.implement_more_comments(" output file for furteher model restart");
  ContinuationFormat.implement_more_comments(" CLX is original 2D restart format, NETCDF extends to 3D simulation");
  ContinuationInterval.implement_more_comments(" Time interval for rolling continuation file updates");
  OutputPath.implement_more_comments(" pathname of output directory\n execution directory is the default \n specify AUTO to automatically create tagged directory:\n output.YYYY.MM.DD.HH.MM");
  LinearSolver.implement_more_comments   (" Solver to be used in implicit linear systems");
  SubLinearSolver.implement_more_comments(" Solver to be used in implicit linear systems, sub-cycling specific (use AUTO except in parallel run)");

  TwoDMometumEqForm.implement_more_comments(" Sould be completed");
  ThreeDPressure.implement_more_comments("test option; not to be used");
  EchoFile.implement_more_comments(" T-UGO setup echoing file");
  MeshFile.implement_more_comments(" Sould be completed");
}

physics_C::physics_C()
{
  gravity_mode.implement ("gravity_mode","<default> CONSTANT","<default> CONSTANT","gravitational constant definition",2,"LATITUDE-DEPENDANT CONSTANT","NONE");
  anelastic_body_tide.implement("anelastic_body_tide","<default> FALSE","<default> FALSE","anelastic body tide flag",2,"FALSE TRUE","NONE");
  gravity_adjustment.implement ("gravity_adjustment","<default> NONE","<default> NONE","gravitational constant adjustment",3,"NONE CONSTANT FORMULA","NONE");
  coriolis_mode.implement("coriolis_mode","<default> LATITUDE-DEPENDANT","<default> LATITUDE-DEPENDANT","coriolis parameter definition",3,"LATITUDE-DEPENDANT CONSTANT BETA-PLAN","NONE");
  coriolis_latitude.implement("coriolis_latitude","<default> 0.0","<default> 0.0","coriolis latitude",0,"NONE","degrees");
  water_viscosity.implement("water_viscosity","<default> 1.0e-06","<default> 1.0e-06","water viscosity",0,"NONE","m²/s");
}

dissipation_C::dissipation_C()
{
   using namespace Keywords ;

   BKGPolygons.implement("BKG_polygones","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   BKGValues.implement("BKG_values","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   BackgroundVelocityFile.implement("background_velocity_file","<default> NONE","<default> NONE","some comments",0,"NONE","File");
   IceShelvesFile.implement("ice_shelves_file","<default> NONE","<default> NONE","ice shelves polygons",0,"FILE","File");
   IceShelvesCoefficient.implement("ice_shelves_coefficients","<default> 2.0","<default> 2.0","friction factor applied below ice shelves",0,"NONE","dimensionless");


   // Make the menu which allows to choose the type of rugosity.
   const std::string initial = KEY_DEFAULT + KEY_QUADRATIC ;
   BottomFrictionType.implement("bottom_friction_type", initial.c_str(), initial.c_str(), "some comments", 3,(KEY_LINEAR + " " + KEY_QUADRATIC + " " + KEY_KARMAN).c_str(),"NONE") ;

   LinearFrictionCoeff.implement("linear_friction_coeff","<default> 2.5e-3","<default> 2.5e-3","default linear bottom friction coefficient (r)",0,"NONE","dimensionless");
   
   QuadraticFrictionCoeff.implement("quadratic_friction_coeff","<default> 2.5e-3","<default> 2.5e-3","default quadratic bottom friction coefficient (Cd)",0,"NONE","dimensionless");
   QFCPolygons.implement("QFC_polygones","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   QFCValues.implement("QFC_values","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   BottomFrictionCoeffFile.implement("bottom_coeff_file","<default> NONE","<default> NONE","prescribed Cd file, LGP1 discretisation",0,"FILE","File");

   KarmanRugosity.implement("karman_rugosity","<default> 1.0e-4","<default> 1.0e-4","some comments",0,"NONE","dimensionless");
   KarmanRugosity.deprecated=true; 

   RoughnessLength.implement("roughness_length","<default> 1.0e-4","<default> 1.0e-4","equivalent roughness length (z0) in quadratic bottom friction coefficient (Cd) computation",0,"NONE","m");
   MinQuadraticFrictionCoeff.implement("min_quadratic_friction_coeff","<default> -999","<default> -999","minimum value for Cd when using roughness formulation",0,"NONE","dimensionless");
   QFCMPolygons.implement("QFCM_polygones","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   QFCMValues.implement("QFCM_values","<default> NONE","<default> NONE","some comments",0,"FILE","File");

   BottomRugosityFile.implement("bottom_rugosity_file","<default> NONE","<default> NONE","prescribed z0 file, LGP1 discretisation",0,"FILE","File");
   RugosityPolygons.implement("rugosity_polygones","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   RugosityValues.implement("rugosity_values","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   RugosityValuesbyRegions.implement("rugosity_values_by_regions","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   
   InternDragAlgorithm.implement("internal_drag_algorithm","<default> 0","<default> 0","some comments",6,"0 1 2 3 4 5","dimensionless");
   InternDragSlope.implement("internal_drag_slope","<default> 0.00","<default> 0.00","some comments",0,"NONE","dimensionless");
   InternDragHmin.implement("internal_drag_hmin","<default> 300","<default> 300","some comments",0,"NONE","dimensionless");
   InternDragHmax.implement("internal_drag_hmax","<default> 500","<default> 500","some comments",0,"NONE","dimensionless");
   IWD01Polygons.implement("internal_drag_polygones","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   IWD01Values.implement("internal_drag_values","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   IWD01ValuesByRegions.implement("internal_drag_values_by_regions","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   InternDragRugosity.implement("internal_drag_rugosity","<default> 0.00","<default> 0.00","some comments",0,"NONE","dimensionless");
   
   MixedLayerCoeff.implement("mixed_layer_coeff","<default> 0.00","<default> 0.00","some comments",0,"NONE","m/s");
   MinBackGroundSpd.implement("min_background_speed","<default> 0.05","<default> 0.05","some comments",0,"NONE","m/s");
   MinhorizontalVisc.implement("min_horizontal_viscosity","<default> 1.00","<default> 1.00","some comments",0,"NONE","m^2/s");
   SmagorinskyCoefficient.implement("smagorinsky_coefficient","<default> 0.28","<default> 0.28","some comments",0,"NONE","dimensionless");
   HorizontalViscMode.implement("horizontal_viscosity_mode","<default> SMAGORINSKY","<default> SMAGORINSKY","some comments",5,"NONE SMAGORINSKY PROPORTIONAL CONSTANT SMAGORINSKY-2D","NONE");
   ShearDragFlag.implement("shear_drag_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");

   BruntVassalaValue.implement("brunt_vassalla_value","<default> 0.002","<default> 0.002","some comments",0,"NONE","s");
   BruntVassalaFile.implement_more_comments("Brûnt-Vasala uniform/default value");

   BruntVassalaFile.implement("brunt_vassalla_file","<default> NONE","<default> NONE","some comments",0,"FILE","File");
   BruntVassalaFile.implement_more_comments("filename (s2r format) with local Brûnt-Vasala values");

   InternDragAlgorithm.implement_more_comments("0: constant tidal excursion\n 1: constant tidal excursion\n 2: wave-dependent tidal excursion\n 3: (1) + latitude  dependant factor" );
} 

tides_C::tides_C()
{
  OnOffFlag.implement("tide_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  AstronomicPotentialFlag.implement("AstronomicPotentialFlag","<default> FALSE","<default> FALSE","enable/disable gravitational forcing",2,"FALSE TRUE","NONE");
  LSAFlag.implement("LSA_flag","<default> FALSE","<default> FALSE","enable/disable LSA forcing",2,"FALSE TRUE","NONE");
  LSADir.implement("LSA_directory","<default> NONE","<default> NONE","pathname for LSA input files",1,"DIRECTORY","Directory");
  LSAConvention.implement("LSA_convention","<default> WAVE.load.nc","<default> WAVE.load.nc","LSA filename formatting",0,"NONE","NONE");
  AtlasDir.implement("Atlas_directory","<default> NONE","<default> NONE","Tidal atlas pathname for online OBC computation",1,"DIRECTORY","Directory");
  AtlasConvention.implement("Atlas_convention","<default> WAVE.tide.nc","<default> WAVE.tide.nc","Atlas filename formatting",0,"NONE","NONE");
  PressureFlag.implement("pressure_flag","<default> FALSE","<default> FALSE","enable/disable atmospheric pressure forcing",2,"FALSE TRUE","NONE");
  PressureDir.implement("pressure_directory","<default> NONE","<default> NONE","pathname for atmospheric pressure input files",1,"DIRECTORY","Directory");
  PressureConvention.implement("pressure_convention","<default> WAVE.tide.nc","<default> WAVE.tide.nc","pressure forcing filename formatting",0,"NONE","NONE");
  admittance.implement("admittance","<default> FALSE","<default> FALSE","enable/disable admittance spectrum extension",2,"FALSE TRUE","NONE");
  equilibrium.implement("equilibrium","<default> FALSE","<default> FALSE","enable/disbale equilibrium tide onnline OBC computation",2,"FALSE TRUE","NONE");
  NodalCorrection.implement("nodal_correction","<default> FALSE","<default> FALSE","enabke/disbale nodal corrections",2,"FALSE TRUE","NONE");
  LSACoeff.implement("LSA_coeff","<default> 0.0","<default> 0.0","some comments",0,"NONE","dimensionless");
  DeformationLoveNum.implement("Deformation_love_num","<default> 0.1","<default> 0.1","linear LSA coefficient",0,"NONE","dimensionless");
  BoundTideFlag.implement("tidal_BC","<default> FALSE","<default> FALSE","enable/disable tidal boundary conditions",2,"TRUE FALSE","NONE");
  BoundTideFile.implement("BC_tide_file","<default> NONE","<default> NONE","tidal boundary conditions input file",1,"FILE","File");
  Spectral2DRun.implement("spectral_solver","<default> FALSE","<default> FALSE","enable/disable 2D-spectral tidal solver",2,"FALSE TRUE","NONE");
  SpectralSolver.implement("solver_type","<default> SpDOMESTIC","<default> SpDOMESTIC","spectral linear solver",6,"AUTO UMFPACK MUMPS PASTIX HIPS SpDOMESTIC","NONE");

  DominantWave_1.implement("dominant_wave_1","<default> M2","<default> M2","1st dominant wave for spectral solver",6,"M2 S2 N2 K1 O1 NONE","NONE");
  DominantWave_2.implement("dominant_wave_2","<default> K1","<default> K1","2nd dominant wave for spectral solver",6,"M2 S2 N2 K1 O1 NONE","NONE");
  PriorSolution.implement("prior_solution","<default> NONE","<default> NONE","prior solution for friction",1,"DIRECTORY","Directory");

  SpectralPaire.implement("spectral_paire","<default> DGP1xLGP2","<default> DGP1xLGP2","discretisation for spectral solver", 9,\
                          "DGP1xLGP2 NCP1xLGP2 DNP1xLGP2 LGP0xLGP1 NCP1xLGP0 NCP1xLGP1 LGP2xLGP2 CQP1xCQP0 CQN1xCQP0","NONE");
  Spectral2DMaxIteration.implement("spectral_max_iteration","<default> 10","<default> 10","max iterations for dominant waves in spectral solver",0,"NONE","integer");
  CompoundMaxIteration.implement("compound_max_iteration","<default> 10","<default> 10","max iterations for compound waves in spectral solver",0,"NONE","integer");
  SpecificData.implement("spectral_specific_data","<default> FALSE","<default> FALSE","use discretisation-specific database for spectral solver",2,"FALSE TRUE","NONE");
  TopoFile.implement("spectral_topo_file","<default> NONE","<default> NONE","specific topography for spectral solver (u-discretisation compatible)",1,"FILE","File");
  SlopeFile.implement("spectral_slope_file","<default> NONE","<default> NONE","specific topography's slope for spectral solver (u-discretisation compatible)",1,"FILE","File");
  SpectralRecycle.implement("spectral_recycle","<default> FALSE","<default> FALSE","recycle spectral solutions for tidal OBCs",2,"FALSE TRUE","NONE");
  ReducedArchive.implement("spectral_reduced_archives","<default> FALSE","<default> FALSE","reduced spectral solutions archive",2,"FALSE TRUE","NONE");
  SpectralArchiveLevel.implement("spectral_archive_level","<default> STANDARD","<default> STANDARD","spectral solutions archive content",4,"MINIMAL STANDARD DEBUGGING EXPLORATION","NONE");
  Spectral3DRun.implement("spectral3D_solver","<default> DISABLED","<default> DISABLED","enable/disable 3D-spectral tidal solver",4,"DISABLED BAROTROPIC BAROCLINIC SW2D+1D","NONE");
  Spectral3DOBCs.implement("spectral3D_OBCs","<default> BAROTROPIC","<default> BAROTROPIC","OBCs",3,"BAROTROPIC BAROCLINIC RADIATIONAL","NONE");
  VDiffusionDuration.implement("diffusion_duration","<default> 72","<default> 72","integration duration turbulence model",0,"NONE","hours");
  VDiffusionTimeStep.implement("diffusion_time_step","<default> 120","<default> 120","time step for turbulence model",0,"NONE","seconds");
  Spectral3DFirstIteration.implement("spectral3D_first_iteration","<default> 0","<default> 0","first iteration for dominant waves in 3D-spectral solver",0,"NONE","integer");
  Spectral3DMaxIteration.implement("spectral3D_max_iteration","<default> 10","<default> 10","max iterations for dominant waves in 3D-spectral solver",0,"NONE","integer");
  active=false;

  ReducedArchive.deprecated=true;

  LSAFlag.implement_more_comments("trigger use of loading/self-attraction forcing in simulation \n" );
  LSADir.implement_more_comments("Directory path for loading/self-attraction input files\n" );
  LSAConvention.implement_more_comments("filename formatting (WAVE will be replaced by wavename (capital spelling), idem for wave(normal spelling)\n");

  AtlasConvention.implement_more_comments("filename formatting convention :\n<WAVE> will automatically be replaced by the adequate wave names (using upper-case spelling).\n idem for <wave> (using lower-case spelling)\n");
  admittance.implement_more_comments("allow admittance technique in forcing and analysis functions \n" );
  equilibrium.implement_more_comments("allow equilibrium value for boundary conditions (long period tides)\n" );

  Spectral2DMaxIteration.implement_more_comments("negative value means dominant waves to be taken from existing archive, frame=-max iterations\n");
  Spectral3DMaxIteration.implement_more_comments("negative value means dominant waves to be taken from existing archive, frame=-max iterations\n");
}

atmosphere_C::atmosphere_C()
{
  OnOffFlag.implement("atmosphere_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  PressureFlag.implement("pressure_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  PressureDir.implement("pressure_directory","<default> AUTO","<default> AUTO","pathname for atmospheric pressure input files",1,"DIRECTORY","Directory");
  PressureConvention.implement("pressure_convention","<default> AUTO","<default> AUTO","pressure forcing filename formatting",0,"NONE","NONE");
  WindFlag.implement("wind_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  LongWave.implement("long_wave_radiation","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ShortWave.implement("short_wave_radiation","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  AtmosDirectory.implement("atmospheric_directory","<default> NONE","<default> NONE","some comments",1,"DIRECTORY","Directory");
  FileFormat.implement("file_format","<default> NETCDF","<default> NETCDF","some comments",5,"NETCDF GRIB CLX ACADEMIC HARMONIC","NONE");

  FileEndian.implement("file_endian","<default> UNIX","<default> UNIX","some comments",2,"UNIX LINUX","NONE");
  FileConvention.implement("file_convention","<default> NONE","<default> NONE","some comments",0,"NONE","NONE");
  FileConvention.implement_more_comments("filename construction rules; example :\n ecmwf-YYYY-MM-DD.nc\n for files named such as ecmwf-2002-01-23.nc ");

  MaskFile.implement("MaskFile","<default> NONE","<default> NONE","some comments",3,"FILE NONE AUTO","File");
  ApplyMask.implement("apply_landsea_mask","<default> FALSE","<default> FALSE","apply land/sea mask on meteorological fields",2,"FALSE TRUE","NONE");
  KeepMeanPressureFlag.implement("keep_mean_pressure","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  KeepMeanWindFlag.implement("keep_mean_wind","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  TimeInterp.implement("time_interpolation","<default> LINEAR","<default> LINEAR","some comments",2,"LINEAR QUADRATIC","NONE");
  BulkFormula.implement("bulk_formula","<default> CLASSIC","<default> CLASSIC","some comments",4,"CLASSIC BULK BULK-REDUCED FRICTION-VELOCITY","NONE");
  WindStress.implement("wind_stress","<default> BULK","<default> BULK","some comments",2,"BULK MODEL","NONE");
  BoundInversebarometer.implement("inversebarometer","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  active=false;
}

OceanWaves_C::OceanWaves_C()
{
  OnOffFlag.implement("OceanWaves_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  Parameterisation.implement("parameterisation","<default> MAC-WILLIAMS","<default> MAC-WILLIAMS","some comments",2,"MAC-WILLIAMS LONGUET-HIGGINS","NONE");
  DragFlag.implement("drag_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  FrictionFlag.implement("friction_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  TimeInterp.implement("time_interpolation","<default> LINEAR","<default> LINEAR","some comments",2,"LINEAR QUADRATIC","NONE");
  Discretisation.implement("discretisation","<default> UNSTRUCTURED","<default> UNSTRUCTURED","some comments",2,"UNSTRUCTURED STRUCTURED","NONE");
  FileDirectory.implement("file_directory","<default> NONE","<default> NONE","some comments",1,"DIRECTORY","Directory");
  FileFormat.implement("file_format","<default> WW3_SG_NETCDF","<default> WW3_SG_NETCDF","input format of ocean waves forcing",4,"WW3_SG_NETCDF WW3_UG_NETCDF WW3_SG_NATIVE WW3_UG_NATIVE","NONE");
  FileConvention.implement("file_convention","<default> NONE","<default> NONE","some comments",0,"NONE","NONE");
  FileConvention.implement_more_comments("filename construction rules; example :\n ecmwf-YYYY-MM-DD.nc\n for files named such as ecmwf-2002-01-23.nc ");

  MaskFile.implement("MaskFile","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  active=false;
}

boundaries_C::boundaries_C()
  {
  type.implement("type","<default> 1","<default> 1","some comments",9,"1 2 3 4 5 7=(1->2) 8 9 10","NONE");
  mode.implement("mode","<default> LOCAL","<default> LOCAL","some comments",2,"LOCAL VARIATIONAL","NONE");
  RelaxTime.implement("relax_time","<default> 2.0","<default> 2.0","some comments",0,"NONE","hours");
  
  ArchivedBCUse.implement("archived_BC_use","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ArchivedBCDirectory.implement("archived_BC_directory","<default> NONE","<default> NONE","some comments",1,"DIRECTORY","directory");
  ArchivedBCFormat.implement("archived_BC_format","<default> BINARY","<default> BINARY","some comments",2,"BINARY NETCDF","NONE");
  ArchivedBCRoot.implement("archived_BC_root","<default> analysis","<default> analysis","some comments",0,"NONE","NONE");
  ArchivedBCUnits.implement("archived_BC_units","<default> MKS","<default> MKS","some comments",2,"MKS CGS","NONE");
  ArchivedBCTransport.implement("archived_BC_transport","<default> u","<default> u","some comments",3,"u U null","NONE");
  
  PeriodicFile.implement("periodic_conditions_file","<default> NONE","<default> NONE","contains instructions for boundary condition periodicity",1,"FILE","File");
  SpongeFlag.implement("sponge_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  PMLfactor.implement("PML_factor","<default> 0.0","<default> 0.0","attenuation factor in PML absorbing zone",0,"NONE","s-1");
  BufferZoneFile.implement("buffer_zone_file","<default> NONE","<default> NONE","contains instructions for open boundaries buffers",1,"FILE","File");
  
  RiversFlag.implement("rivers_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  RiversFile.implement("rivers_file","<default> rivers.dat","<default> rivers.dat","some comments",1,"FILE","File");
  RiversTimeTemplate.implement("rivers_time_template","<default> CNESTIME","<default> CNESTIME","some comments",0,"NONE","NONE");
  SolidCondition.implement("solid_boundary_condition","<default> FREE-SLIP","<default> FREE-SLIP","some comments",2,"FREE-SLIP CONDITIONAL-NO-SLIP","NONE");
  NoSlipSize.implement("no_slip_size","<default> 100.0","<default> 100.0","some comments",0,"NONE","m");

  ArchivedBCUse.implement_more_comments("Open boundary conditions taken from a specified large scale model (nesting)\n"
                                        "nesting available for T-UGOm type of archives" );

  ArchivedBCDirectory.implement_more_comments("Path to the directory containing the specified large scale model archives\n");

  ArchivedBCRoot.implement_more_comments("rootname for the specified large scale model archive files"
                                         "model will seek for ROOTNAME.YYYY.MM" );

  ArchivedBCUnits.implement_more_comments("units system of the specified large scale model archive files"
                                          "here to keep compatibility with MOG2D archives (CGS)" );

  ArchivedBCTransport.implement_more_comments("velocity at boundaries : use zero, extract u or transport U\n"
                                              "extraction taken from the specified large scale model \n"
                                              "--------------------------------------- \n"
                                              "needed in characteristics boundary conditions (Flather)" );

  RiversFlag.implement_more_comments("trigger river discharge in simulation\n");
  RiversFile.implement_more_comments("contains river discharge informations\n");
  
  RiversTimeTemplate.deprecated=true;

  }

ice_C::ice_C()
{
  OnOffFlag.implement("ice_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  IceCoverFile.implement("ice_cover_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  IceElasticityFile.implement("ice_elasticity_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  TimeStep.implement("TimeStep","<default> 60","<default> 60","some comments",0,"NONE","s");
  MomentumSolver.implement("MomentumSolver","<default> 1","<default> 1","some comments",4,"0 1 2 3","NONE");
  Rheology.implement("Rheology","<default> 1","<default> 1","some comments",2,"0 1","NONE");
  active=false;
}

topography_C::topography_C() {
  OnOffFlag.implement("OnOffFlag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  BottomMotionFlag.implement("bottom_motion_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  BottomMotionDirectory.implement("bottom_motion_directory","<default> NONE","<default> NONE","some comments",1,"DIRECTORY","Directory");
  BottomMotionFileConvent.implement("bottom_motion_file_convention","<default> NONE","<default> NONE","some comments",0,"NONE","NONE");

  TopoFile.implement("topo_file","<default> NONE","<default> NONE","bathymetry file",1,"FILE","File");
  SlopeFile.implement("slope_file","<default> NONE","<default> NONE","bathymetry gradient's file",1,"FILE","File");
  ZminFile.implement("zmin_file","<default> NONE","<default> NONE","minimum elevation file",1,"FILE","File");
  RelativeMinDepth.implement("relative_minimum_depth","<default> 5.0","<default> 5.0","minimum depth relative to minimum elevation",0,"NONE","metres");
  MinDepth.implement("minimum_depth","<default> 10.0","<default> 10.0","overall minimum depth in model bathymetry",0,"NONE","metres");
  DryingMinDepth.implement("drying_minimum_depth","<default> 0.10","<default> 0.10","numerical minimum depth set at dry nodes",0,"NONE","metres");

  DepthOffsetPolygons.implement("depth_offset_polygons","<default> NONE","<default> NONE","depth offset file",1,"FILE","File");
  DepthOffsetValues.implement("depth_offset_values","<default> NONE","<default> NONE","depth offset file",1,"FILE","File");
  MeanLevelFile.implement("mean_level_file","<default> NONE","<default> NONE","depth offset file",1,"FILE","File");
  DepthScaleFactor.implement("depth_scale_factor","<default> 1.0","<default> 1.0","depth scale factor",0,"NONE","dimensionless");
 
  active=false;
}

misc_C::misc_C() {
  NonTidalLoadFlag.implement("non_tidal_loading_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  SmoothingFlag.implement("smoothing_flag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  checks.implement("checks","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  Check_Wconsistency.implement("Check_Wconsistency","<default> FALSE","<default> FALSE","some comments",2,"TRUE FALSE","NONE");
  CompatibilityScale.implement("CompatibilityScale","<default> 1.0","<default> 1.0","some comments",0,"NONE","dimensionless");

  BarotropicW.implement("barotropic_w","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  BarotropicW.deprecated=true;

  coupling_2D_3D.implement("2D_3D_coupling","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");

  StabilityCheck.implement("Stability_Check","<default> TRUE","<default> TRUE","some comments",2,"TRUE FALSE","NONE");
  StabilityControl.implement("Stability_Control","<default> DYNAMIC","<default> DYNAMIC","some comments",2,"DYNAMIC STATIC","NONE");
  SubcyclingCFLratio.implement("SubcyclingCFLratio","<default> 2","<default> 2","some comments",0,"NONE","NONE");

  MaxDeltaElevation.implement("Max_Delta_Elevation","<default> 4.0","<default> 4.0","some comments",0,"NONE","m");
  MaxDeltaCurrents.implement("Max_Delta_Currents","<default> 4.0","<default> 4.0","some comments",0,"NONE","m/s");
  MaxDeltaPeriod.implement("Max_Delta_Period","<default> 3600.","<default> 3600.","some comments",0,"NONE","s");
  nSubCycles.implement("nSubCycles","<default> 1","<default> 1","some comments",0,"NONE","NONE");
  SubcyclingPersistence.implement("persistence","<default> 6.","<default> 6.","some comments",0,"NONE","hours");
  SubcyclingTemporization.implement("SubcyclingTemporization","<default> 10","<default> 10","some comments",0,"NONE","NONE");
  SubcyclingExtension.implement("SubcyclingExtension","<default> 5","<default> 5","some comments",0,"NONE","NONE");
  SubcyclingLockLimit.implement("SubcyclingLockLimit","<default> 5","<default> 5","some comments",0,"NONE","NONE");
  SmoothZField.implement("SmoothZField","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  SmoothUField.implement("SmoothUField","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  SmoothAdvection.implement("SmoothAdvection","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  SmoothDiffusion.implement("SmoothDiffusion","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  FlushTracking.implement("flush_tracking","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ToponymsFile.implement("toponyms_file","<default> NONE","<default> NONE","ocean partition codes",0,"FILE","File");
}

fearchive_C::fearchive_C() {
  OnOffFlag.implement("OnOffFlag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  StateFile.implement("state_file","<default> tugo.state","<default> tugo.state","some comments",1,"FILE","File");
  VectorFile.implement("vector_file","<default> tugo.forcing","<default> tugo.forcing","some comments",1,"FILE","File");
  ArchFormat.implement("archive_format","<default> NETCDF","<default> NETCDF","some comments",2,"NETCDF CLX","NONE"); 
  ArchConvention.implement("archive_convention","<default> NONE","<default> NONE","some comments",0,"NONE","NONE");
  ArchStart.implement("archive_start","<default> 0.0","<default> 0.0","some comments",0,"NONE","days");
  ArchInterval.implement("archive_interval","<default> 10800.0","<default> 10800.0","some comments",0,"NONE","seconds");
  ArchBits.implement("archive_bits","<default> 2","<default> 2","some comments",3,"1 2 4","NONE"); 
  Arch_huv_Scale.implement("archive_huvscale","<default> 100 100 100","<default> 100 100 100","some comments",0,"NONE","factor");
  Arch_pW_Scale.implement("archive_pWsxWsy_scale","<default> 1000 1000 1000","<default> 1000 1000 1000","some comments",0,"NONE","factor");

  ArchOverAppendFlag.implement("archive_overwright_append_flag","<default> APPEND","<default> APPEND","some comments",2,"APPEND OVERWRITE","NONE");
  ArchOverAppendFlag.deprecated=true;

  OverwriteAppendFlag.implement("overwrite_append_flag","<default> APPEND","<default> APPEND","some comments",2,"APPEND OVERWRITE","NONE");

  z_flag.implement("surface_elevation","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  ubar_flag.implement("barotropic_currents","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  horizontal_diffusion2D.implement("horizontal_diffusion2D","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  horizontal_advection2D.implement("horizontal_advection2D","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  transport2D.implement("transport2D","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  divergence2D.implement("divergence2D","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  u_flag.implement("horizontal_currents","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  w_flag.implement("vertical_velocity","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  T_flag.implement("temperature","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  S_flag.implement("salinity","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  R_flag.implement("density","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  horizontal_diffusion.implement("horizontal_diffusion","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  horizontal_advection.implement("horizontal_advection","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  Cd_flag.implement("friction_coefficient","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ws_flag.implement("wind_stress","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  Pa_flag.implement("atmospheric_surface_pressure","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ibd_flag.implement("IB-departure","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  lateral_mixing.implement("lateral_mixing","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  radiation_stress.implement("radiation_stress","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  radiation_stress.implement_more_comments("activate writing of radiation stress components (Sxx, Sxy, Syy) from Wave model");
  gradSxx.implement("gradSxx","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  gradSxx.implement_more_comments("activate writing of gradSxx");
  active=false;
}

structured_C::structured_C() {
  OnOffFlag.implement("OnOffFlag","<default> 0","<default> 0","some comments",3,"0 1 2","NONE");
  zone.implement("zone","<default> AUTO","<default> AUTO","some comments",3,"AUTO FILE SYMPHO","NONE");
  OverwriteAppendFlag.implement("overwrite_append_flag","<default> APPEND","<default> APPEND","some comments",2,"APPEND OVERWRITE","NONE");
  Start.implement("start","<default> 0.0","<default> 0.0","some comments",0,"NONE","days");
  Interval.implement("interval","<default> 10800.0","<default> 10800.0","some comments",0,"NONE","seconds");
  ArchiveStart.implement("archive_start","<default> 0.0","<default> 0.0","some comments",0,"NONE","days");
  ArchiveStart.deprecated=true;
  ArchiveInterval.implement("archive_interval","<default> 10800.0","<default> 10800.0","some comments",0,"NONE","seconds");
  ArchiveInterval.deprecated=true;
  z_flag.implement("surface_elevation","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  zrange.implement("surface_elevation_extrema","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  MKE2D.implement("mean_kinetic_energy2D","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ubar_flag.implement("barotropic_currents","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  horizontal_diffusion2D.implement("horizontal_diffusion2D","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  horizontal_advection2D.implement("horizontal_advection2D","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  transport2D.implement("transport2D","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  divergence2D.implement("divergence2D","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  Cd_flag.implement("friction_coefficient","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  z0_flag.implement("bottom_rugosity","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  Pa_flag.implement("atmospheric_surface_pressure","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ws_flag.implement("wind_stress","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ibd_flag.implement("IB-departure","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  u_flag.implement("horizontal_currents","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  w_flag.implement("vertical_velocity","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  T_flag.implement("temperature","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  S_flag.implement("salinity","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  R_flag.implement("density","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  horizontal_diffusion.implement("horizontal_diffusion","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  horizontal_advection.implement("horizontal_advection","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  pressure_gradient.implement("pressure_gradient","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  lateral_mixing.implement("lateral_mixing","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  tracer2D.implement("tracer_2d","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  radiation_stress.implement("radiation_stress","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  radiation_stress.implement_more_comments("activate writing of radiation stress components (Sxx, Sxy, Syy) from Wave model");
}

sample_C::sample_C() {
  OnOffFlag.implement("OnOffFlag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  file.implement("file","<default> tugoPlot.input","<default> tugoPlot.input","some comments",1,"FILE","File");
  start.implement("start","<default> 0.0","<default> 0.0","some comments",0,"NONE","seconds");
  SaveInterval.implement("save_interval","<default> 1800","<default> 1800","some comments",0,"NONE","seconds");
  OverAppendFlag.implement("overwright_append_flag","<default> APPEND","<default> APPEND","some comments",2,"APPEND OVERWRITE","NONE");
  active=false;
  SectionOnOffFlag.implement("SectionOnOffFlag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  SectionInputFile.implement("SectionInputFile","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  SectionStart.implement("SectionStart","<default> 0.0","<default> 0.0","some comments",0,"NONE","seconds");
  SectionSaveInterval.implement("SectionSaveInterval","<default> 1800","<default> 1800","some comments",0,"NONE","seconds");
  SectionOverAppendFlag.implement("SectionOverAppendFlag","<default> APPEND","<default> APPEND","some comments",2,"APPEND OVERWRITE","NONE");
  SectionActive=false;
}

analysis_C::analysis_C() {
  OnOffFlag.implement("OnOffFlag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ListFile.implement("list_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  AutoCompletion.implement("auto_completion","<default> TRUE","<default> TRUE","if list file not given, do spectrum auto-completion",2,"FALSE TRUE","NONE");
  HarmoStart.implement("harmonic_start","<default> 2*#MODEL.$spinup./3600/24.0","<default> 2*#MODEL.$spinup./3600/24.0","some comments",0,"NONE","days");
  SampleInterval.implement("sample_interval","<default> 1800","<default> 1800","some comments",0,"NONE","seconds");
  ComputeInterval.implement("compute_interval","<default> 365","<default> 365","some comments",0,"NONE","days");
  Strategy.implement("strategy","<default> FULL","<default> FULL","some comments",3,"FULL REINITIALISE SEQUENTIAL","NONE");
  SequentialRecycle.implement("sequential_recycle","<default> MANAGED","<default> MANAGED","recycle sequential solutions for tidal OBCs",3,"MANAGED FALSE TRUE","NONE");
  SequentialArchiveLevel.implement("sequential_archive_level","<default> STANDARD","<default> STANDARD","sequential solutions archive content",2,"MINIMAL STANDARD","NONE");
  LGP1ArchiveFlag.implement("LGP1ArchiveFlag","<default> TRUE","<default> TRUE","create a LGP1xLGP1 ASCII archive for compatibility purposes",2,"TRUE FALSE","NONE");
  NETCDFArchiveFlag.implement("NETCDFArchiveFlag","<default> TRUE","<default> TRUE","create a netcdf UG archive (on computational discretisation)",2,"TRUE FALSE","NONE");
  grid.implement("zone","<default> AUTO","<default> AUTO","some comments",3,"AUTO FILE SYMPHO","NONE");
  SGArchiveFlag.implement("SGFlag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  z_flag.implement("surface_elevation","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  u_flag.implement("horizontal_currents","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  LSA_flag.implement("LSA_potential","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  active=false;
  SGsave=false;
  ListFile.implement_more_comments("contains additional information on tidal spectrum (forcing and analysis) \n" );
}

energy_C::energy_C() {
  OnOffFlag.implement("OnOffFlag","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  start.implement("start","<default> 0.0","<default> 0.0","some comments",0,"NONE","units to be verificated");
  interval.implement("interval","<default> 1800","<default> 1800","some comments",0,"NONE","units to be verificated");
  TimeStepComp.implement("time_step_energy_budget","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  HarmoComp.implement("harmonic_energy_budget","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  HarmoCompByRegions.implement("harmonic_energy_budget_by_regions","<default> NONE","<default> NONE","some comments",0,"FILE","File");
  ResonanceExploration.implement("resonance_exploration","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ResonanceFrequencies.implement("resonance_frequencies","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  active=false;
}

tracers2D_C::tracers2D_C() {
  OnOffFlag.implement("tracers computation","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  start.implement("start","<default> 0.0","<default> 0.0","some comments",0,"NONE","units to be verificated");
  interval.implement("interval","<default> 1800","<default> 1800","some comments",0,"NONE","seconds");
  TimeIntegration.implement("time_integration","<default> LEAPFROG","<default> LEAPFROG","some comments",3,"FORWARD LEAPFROG EULER-BF","NONE");
  InputFile.implement("input_file","<default> tracers.dat","<default> tracers.dat","some comments",1,"FILE","File");
  InitialisationFile.implement("initialisation_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  HorizontalDiffusionMode.implement("horizontal_diffusion_mode", "<default> DEFAULT","<default> DEFAULT","some comments",5,"DEFAULT NONE SMAGORINSKY CONSTANT UPWINDING","NONE");
  HorizontalDiffusionCoefficient.implement("horizontal_diffusion_coefficient","<default> 0.0","<default> 0.0","some comments",0,"NONE","m^2/s");
  SmagorinskyCoefficient.implement("smagorinsky_coefficient","<default> 2.8","<default> 2.8","some comments",0,"NONE","dimensionless");
  UpwindCoefficient.implement("upwind_coefficient","<default> 0.1","<default> 0.1","some comments",0,"NONE","dimensionless");
  active=false;

  HorizontalDiffusionMode.implement_more_comments(" Controls 2D tracers diffusion\n"
                                                  " DEFAULT means identical to the momentum diffusion \n"
                                                  " UPWIND is diffusion tuned for advection upwinding \n");
  HorizontalDiffusionCoefficient.implement_more_comments(" 2D tracers diffusion coefficient:\n"
                                                  " value used for constant diffusion coefficient scheme \n"
                                                  " minimum value in the other schemes, except NONE \n");
  UpwindCoefficient.implement_more_comments(" 2D tracers upwinding diffusion coefficient:\n"
                                            " tuning value in upwinding diffusion coefficient scheme \n"
                                            " ranges from 0 to 1 \n");
}

drifters2D_C::drifters2D_C() {
  OnOffFlag.implement("drifters2D computation","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  start.implement("start","<default> 0.0","<default> 0.0","some comments",0,"NONE","units to be verificated");
  interval.implement("interval","<default> 1800","<default> 1800","some comments",0,"NONE","units to be verified");
  InputFile.implement("input_file","<default> drifters2D.dat","<default> drifters2D.dat","some comments",1,"FILE","File");
  active=false;
} 

dynamics2D_C::dynamics2D_C() {
  WaveEquationForm.implement("WaveEquationForm","<default> GWE-CLX","<default> GWE-CLX","continuity equation type",5,"GWE-CLX GWE-XXX WE-EXPLICIT WE-IMPLICIT FV-IMPLICIT","NONE");
  Continuity2DIntegration.implement("Continuity2D_integration","<default> FORWARD","<default> FORWARD","continuity time scheme",3,"FORWARD LEAPFROG UNSPECIFIED","NONE");
  Momentum2DIntegration.implement("Momentum2D_integration","<default> FORWARD","<default> FORWARD","momentum time scheme",3,"FORWARD LEAPFROG EULER-BF","NONE");
  Momentum2DMode.implement("Momentum2D_mode","<default> VELOCITY","<default> VELOCITY","velocity/transport option in model state variables",2,"VELOCITY TRANSPORT","NONE");
  Advection2Dmode.implement("Advection2D_mode","<default> NON-CONSERVATIVE","<default> NON-CONSERVATIVE","some comments",3,"NON-CONSERVATIVE CONSERVATIVE LAGRANGIAN","NONE");
  AdvectionDiffusionDiscretisation.implement("AdvectionDiffusion_discretisation","<default> NODE-WISE","<default> NODE-WISE","some comments",2,"NODE-WISE ELEMENT-WISE","NONE");
  SLVDiscretisation.implement("SLV_discretisation","<default> LGP1-CG","<default> LGP1-CG","elevation/pressure/tracer discretisation",6,"LGP0 LGP1-CG LGP1-DG LGP2-CG CQP0 CQP1","NONE");
  V2DDiscretisation.implement("v2D_discretisation","<default> NCP1-CG","<default> NCP1-CG","velocity discretisation",8,"LGP0 LGP1 LGP1-DG NCP1-CG NCP1-DG CQP0 CQP1 CQN1","NONE");
  G2DDiscretisation.implement("G2D_discretisation","<default> INTRINSIC","<default> INTRINSIC","pressure gradient discretisation",8,"INTRINSIC AUTO LGP0 LGP1 LGP1-DG NCP1-CG NCP1-DG LGP2","NONE");
  WaveScalarProduct.implement("wave_scalar_prouct","<default> QUADRATURE","<default> QUADRATURE","type of scalar product in wave equation discretization",2,"QUADRATURE INTEGRALE","NONE");
  tau0.implement ("tau_0","<default> 2e-4","<default> 2e-4","weight of CE in GWE (more by clicking...)",0,"NONE","s-1");
  theta.implement("theta","<default> 0.75","<default> 0.75","implicit/explicit coefficient (more by clicking...)",0,"NONE","units to be verificated");
  BarotropicW.implement("barotropic_w","<default> FALSE","<default> FALSE","compute vertical velocity checks",2,"FALSE TRUE","NONE");
  HorizontalViscMode.implement("horizontal_viscosity_mode","<default> SMAGORINSKY","<default> SMAGORINSKY","horizontal diffusion ",5,"NONE CONSTANT PROPORTIONAL SMAGORINSKY UPWINDING","NONE");
  UpwindCoefficient.implement("upwind_coefficient","<default> 0.1","<default> 0.1","upwind in diffusion computation(more by clicking...)",0,"NONE","dimensionless");
  AsselinFilteringH.implement("AsselinFilteringH","<default> FALSE","<default> FALSE","Asselin coefficient for elevation (Leapfrog)",2,"TRUE FALSE","NONE");
  AsselinFilteringU.implement("AsselinFilteringU","<default> TRUE","<default> TRUE","Asselin coefficient for velocity (Leapfrog)",2,"TRUE FALSE","NONE");
  AsselinCoefficientC1.implement("AsselinCoefficientC1","<default> 0.1","<default> 0.1","c Xn + 1/2*c(Xn-1 + Xn+1)",0,"NONE","dimensionless");

  tau0.implement_more_comments (" Continuity equation weight in the generalized wave equation");
  theta.implement_more_comments(" When using semi-implicit, time centering coefficient\n 0 is explicit, 1 fully implicit");

  UpwindCoefficient.implement_more_comments(" 2D momentum upwinding diffusion coefficient:\n"
                                            " tuning value in upwinding diffusion coefficient scheme \n"
                                            " ranges from 0 to 1 \n");
  }

dynamics3D_C::dynamics3D_C() {
  OnOffFlag.implement("3D_computation","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  TimeStep.implement("3D_time_step","<default> $2D_time_step","<default> $2D_time_step","some comments",0,"NONE","seconds");
  MomentumIntegration.implement("Momentum_integration","<default> LEAPFROG","<default> LEAPFROG","some comments",3,"FORWARD LEAPFROG EULER-BF","NONE");
  TracersIntegration.implement("Tracers_integration","<default> FORWARD","<default> FORWARD","some comments",3,"FORWARD LEAPFROG EULER-BF","NONE");
  TracersIntegration.deprecated=true;
  RhoMode.implement("3D_mode","<default> BAROCLINIC","<default> BAROCLINIC","some comments",2,"BAROCLINIC BAROTROPIC","NONE");
  TracersMode.implement("tracers_mode","<default> DIAGNOSTIC","<default> DIAGNOSTIC","some comments",3,"STATIC DIAGNOSTIC PROGNOSTIC","NONE");
  TracersMode.deprecated=true;
  PressureMode.implement("pressure_mode","<default> HYDROSTATIC","<default> HYDROSTATIC","some comments",2,"HYDROSTATIC NON-HYDROSTATIC","NONE");
  WDerivation.implement("w_derivation","<default> DIRECT","<default> DIRECT","some comments",2,"DIRECT WE-CONSISTENT","NONE");
  VDiscretisation.implement("v_discretisation","<default> NCP1","<default> NCP1","some comments",1,"NCP1","NONE");
  WDiscretisation.implement("w_discretisation","<default> QLP1","<default> QLP1","some comments",2,"QLP0 QLP1","NONE");
  TracersDiscretisation.implement("tracers_discretisation","<default> LGP0","<default> LGP0","some comments",3,"LGP0 LGP1 QLP1","NONE");
  TracersDiscretisation.deprecated=true;
  TracersAdvectionCentering.implement("tracers_advection_centering","<default> 0.5","<default> 0.5","some comments",0,"NONE","dimensionless");
  TracersAdvectionCentering.deprecated=true;
  TracerInitialisation.implement("tracer_initialisation","<default> HOMOGENEOUS","<default> HOMOGENEOUS","some comments",3,"HOMOGENEOUS PROFILE MODEL","NONE");
  TracerInitialisation.deprecated=true;
  InitialisationDirectory.implement("initialisation_files_directory","<default> NONE","<default> NONE","some comments",1,"DIRECTORY","Directory");
  InitialisationDirectory.deprecated=true;
  TFile.implement("temperature_initialisation_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  TFile.deprecated=true;
  SFile.implement("salinity_initialisation_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  SFile.deprecated=true;
  RFile.implement("density_initialisation_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  RFile.deprecated=true;
  PFile.implement("tracer_profile_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  PFile.deprecated=true;
  NLayersMax.implement("max_number_layers","<default> 20","<default> 20","some comments",0,"NONE","NONE");
  NUpperHLayers.implement("number_of_upper_horizontal_layers","<default> 0","<default> 0","some comments",0,"NONE","NONE");
  UpperHLayerLimit.implement("limit_of_upper_horizontal_layers","<default> -50","<default> -50","some comments",0,"NONE","m");
  LevelDistribution.implement("level_distribution","<default> UNIFORM","<default> UNIFORM","some comments",4,"COSINUS BOTTOM-COSINUS UNIFORM FILE","NONE");
  LevelDisplacement.implement("level_displacement","<default> SURFACE-ONLY","<default> SURFACE-ONLY","some comments",3,"SURFACE-ONLY ALL RIGID-LID","NONE");
  TurbulenceClosure.implement("turbulence_closure","<default> CONSTANT","<default> CONSTANT","some comments",7,"CONSTANT YAMADA GALPERIN GASPAR TKE KOMEGA HOMOGENEOUS","NONE");
  LateralDiffusion.implement("lateral_diffusion_coefficient","<default> 0.0","<default> 0.0","some comments",0,"NONE","m^2/s");
  VerticalDiffusion.implement("vertical_diffusion_coefficient","<default> 1.e-04","<default> 1.e-04","some comments",0,"NONE","m^2/s");
  TracersVerticalDiffusion.implement("tracer_vertical_diffusion_coefficient","<default> 1.e-04","<default> 1.e-04","some comments",0,"NONE","m^2/s");
  TracersVerticalDiffusion.deprecated=true; 
}

tracers3D_C::tracers3D_C() {
  OnOffFlag.implement("3D_computation","<default> TRUE","<default> TRUE","some comments",2,"FALSE TRUE","NONE");
  TracersIntegration.implement("Tracers_integration","<default> FORWARD","<default> FORWARD","some comments",3,"FORWARD LEAPFROG EULER-BF","NONE");
  TracersMode.implement("tracers_mode","<default> DIAGNOSTIC","<default> DIAGNOSTIC","some comments",3,"STATIC DIAGNOSTIC PROGNOSTIC","NONE");
  TracersDiscretisation.implement("tracers_discretisation","<default> LGP0","<default> LGP0","some comments",3,"LGP0 LGP1 QLP1","NONE");
  TracersAdvectionCentering.implement("tracers_advection_centering","<default> 0.5","<default> 0.5","some comments",0,"NONE","dimensionless");
  TracerInitialisation.implement("tracer_initialisation","<default> HOMOGENEOUS","<default> HOMOGENEOUS","some comments",5,"HOMOGENEOUS UNIFORM-N PROFILE UG-ARCHIVE SG-ARCHIVE","NONE");
  InitialisationDirectory.implement("initialisation_files_directory","<default> NONE","<default> NONE","some comments",1,"DIRECTORY","Directory");
  TFile.implement("temperature_initialisation_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  TVariable.implement("temperature_variable","<default> AUTO","<default> AUTO","name of temperature variable",0,"NONE","dimensionless");
  SFile.implement("salinity_initialisation_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  SVariable.implement("salinity_variable","<default> AUTO","<default> AUTO","name of salinity variable",0,"NONE","dimensionless");
  RFile.implement("density_initialisation_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  RVariable.implement("density_variable","<default> NONE","<default> NONE","name of density variable",0,"NONE","dimensionless");
  PFile.implement("tracer_profile_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  TurbulenceClosure.implement("turbulence_closure","<default> CONSTANT","<default> CONSTANT","some comments",7,"CONSTANT YAMADA GALPERIN GASPAR TKE KOMEGA HOMOGENEOUS","NONE");
  LateralDiffusion.implement("tracers_lateral_diffusion_coefficient","<default> 0.0","<default> 0.0","some comments",0,"NONE","m^2/s");
  HorizontalDiffusion.implement("tracers_horizontal_diffusion_coefficient","<default> 0.0","<default> 0.0","some comments",0,"NONE","m^2/s");
  VerticalDiffusion.implement("tracer_vertical_diffusion_coefficient","<default> 1.e-04","<default> 1.e-04","some comments",0,"NONE","m^2/s");
  HorizontalSmooth.implement("tracers_Hsmoothing","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ComputeDiagnostics.implement("compute_diagnostics","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
}

streamf_C::streamf_C() {
  OnOffFlag.implement("stream_function_computation","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ComputeInterval.implement("compute_interval","<default> 1800","<default> 1800","some comments",0,"NONE","units to be verified");
  active=false;
}

constraint_C::constraint_C() {
  OnOffFlag.implement("stream_function_computation","<default> FALSE","<default> FALSE","some comments",2,"FALSE TRUE","NONE");
  ElevationConstraintFile.implement("elevation_constraint_file","<default> NONE","<default> NONE","some comments",1,"FILE","File");
  active=false;
}

