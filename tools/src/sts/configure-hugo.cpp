
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

#include "tugo-prototypes.h"
#include "dry.h"
#include "ice.h"
#include "poc-netcdf.hpp"
#include "matrix.h"

#include "init_config.hpp"
#include "keywords.hpp"
#include <string>
#include <sys/stat.h> //mkdir
#include "exceptions.hpp"

extern void import_nl_dot_inp_file(tugo_cfg_C *tugo_cfg,char *filename);
extern  int gradient_natural_discretisation(int discretisation);
 
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int create_OuputDirectory(tugo_cfg_C *tugo_cfg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace Keywords;

  int i, j, n = 499, Days, Hours, Minutes, idum, ndum, hour, nitems;
  int status;

  size_t l;
  char msg[1024];
  char key[500];
  char dum[500];
  char *s;
  const char *CommentsEnd = "EndComments";
  float Seconds, fdum1, dummy, fdum;
  struct stat mystat;
  int ErrorNumber;
  
  if(strcmp(gOutputPath,"AUTO")==0) {
    time_t creation_time;
    tm *now;
    free(gOutputPath);
    status = time(&creation_time);
    now=localtime(&creation_time);
    asprintf(&gOutputPath,"./output-%4.4d.%2.2d.%2.2d-%2.2d.%2.2d",now->tm_year+1900,now->tm_mon+1,now->tm_mday,now->tm_hour,now->tm_min);
    status=mkdir(gOutputPath,0777);
    }
  else {
    status=stat(gOutputPath,&mystat);
    if(status==-1) {
      status=mkdir(gOutputPath,0777);
      }
    }
  
  printf(SEPARATOR_1);
  printf("Path for simulation outputs  : %s \n\n", gOutputPath) ;

  status=stat(gOutputPath,&mystat);

  return(status);

//  ErrorNumber = geterrno();
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void complete_configuration01(tugo_cfg_C *tugo_cfg, int create)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace Keywords ;

  int i, j, k, n = 499;
  int Days, Hours, Minutes, idum, ndum, hour, nitems;
  int status;

  size_t l;
  char msg[1024],tmp[1024];
  char key[500];
  char *s, *creation_date;
  const char *CommentsEnd = "EndComments";
  float Seconds, fdum1, dummy, fdum;
  //char *pointer;
  char pointer [100] ;
  time_t creation_time;
  string mystring,setting;

  status = time(&creation_time);

  creation_date=poc_strdup(ctime(&creation_time));

  if(create) status=create_OuputDirectory(tugo_cfg);

  mystring=gOutputPath + (string) "/" + tugo_cfg->model->EchoFile.face_value();
  sscanf(mystring.c_str(), "%s", tmp) ;
  
/**----------------------------------------------------------------------------
  now MPI safe */
  if(gnCPUs==1) {
    strcpy(gRecordFile,tmp);
    }
  else {
    sprintf(gRecordFile,"%s-%4.4d",tmp,gCPU_ID);
    }
    
  if(echo!=NULL) {
    fclose(echo);
    }
    
  echo = fopen(gRecordFile, "w");
  if(echo == NULL) {
    echo = fopen("/dev/null", "w");
//     sprintf(msg, "Cannot open echo file %s \n", gRecordFile);
//     check_error(-1, msg, __LINE__, __FILE__, 1);
    }

  fprintf(echo, "Simulation using %s\n", tugo_version);
  /* STS: NO hg_version */
  fprintf(echo, "Simulation launched %s\n", creation_date);
  fprintf(echo, "Record file used:- %s\n", gRecordFile);
  fprintf(echo, "Path for simulation outputs  : %s \n", gOutputPath) ;

/*----------------------------------------------------------------------
  input files section */
  if(!tugo_cfg->tides->active) {
    tugo_cfg->tides->AstronomicPotentialFlag.actual_value.assign("FALSE");
    tugo_cfg->tides->LSAFlag.actual_value.assign("FALSE");
    tugo_cfg->tides->BoundTideFlag.actual_value.assign("FALSE");
    }

  strcpy(gMeshFileFormat, tugo_cfg->model->MeshFormat.face_value().c_str()) ;
  fprintf(echo, " Geometry file format:-  %s\n", gMeshFileFormat) ;

  if(tugo_cfg->model->MeshType.face_value() == (string) "TRIANGULAR") {
    gMeshFileType=FE_TRIANGLE;
    }
  else if(tugo_cfg->model->MeshType.face_value() == (string) "QUADRANGULAR") {
    gMeshFileType=FE_QUADRANGLE;
    }
  else if(tugo_cfg->model->MeshType.face_value() == (string) "STRUCTURED") {
    gMeshFileType=FE_STRUCTURED;
    }
  else {
    gMeshFileType=FE_UNDEFINED;
    }
  fprintf(echo, " Geometry file type:-  %s\n", tugo_cfg->model->MeshType.face_value().c_str()) ;

  strcpy(GridNameRoot, tugo_cfg->model->rootname.face_value().c_str()) ;
  fprintf(echo, " Root name of grid geometry files:- %s\n", GridNameRoot) ;

  strcpy(GridName, tugo_cfg->model->MeshFile.face_value().c_str()) ;
  fprintf(echo, " Alternative name of grid geometry files:- %s\n", GridName) ;

  strcpy(BelFile, tugo_cfg->model->BelFile.face_value().c_str()) ;
  fprintf(echo, " Boundary (.bel) file:- %s\n", BelFile);

  strcpy(PeriodicFile, tugo_cfg->boundaries->PeriodicFile.face_value().c_str()) ;
  fprintf(echo, " Boundary periodicity file:- %s\n", PeriodicFile);

  strcpy(WaveFile, tugo_cfg->analysis->ListFile.face_value().c_str()) ;
  fprintf(echo, " Optional wave list file:- %s\n", WaveFile);

  strcpy(IceCoverFile, tugo_cfg->ice->IceCoverFile.face_value().c_str()) ;
  fprintf(echo, " Optional ice cover file:- %s\n", IceCoverFile);

  strcpy(IceElasticityFile, tugo_cfg->ice->IceElasticityFile.face_value().c_str()) ;
  fprintf(echo, " Optional ice elasticity file:- %s\n", IceElasticityFile);

  strcpy(BackgroundFile, tugo_cfg->dissipation->BackgroundVelocityFile.face_value().c_str()) ;
  fprintf(echo, " Optional background velocity file:- %s\n", BackgroundFile);

  strcpy(FrictionFile, tugo_cfg->dissipation->BottomFrictionCoeffFile.face_value().c_str()) ;
  fprintf(echo, " Optional friction coefficient file:- %s\n", FrictionFile);

  strcpy(TopoFile, tugo_cfg->topography->TopoFile.face_value().c_str());
  fprintf(echo, " Optional topography file:- %s\n", TopoFile);

  strcpy(SlopeFile, tugo_cfg->topography->SlopeFile.face_value().c_str()) ;
  fprintf(echo, " Optional slope file:- %s\n", SlopeFile);

  strcpy(ZminFile, tugo_cfg->topography->ZminFile.face_value().c_str()) ;
  fprintf(echo, " Optional zmin file:- %s\n", ZminFile);

  strcpy(RugosityFile, tugo_cfg->dissipation->BottomRugosityFile.face_value().c_str()) ;
  fprintf(echo, " Optional rugosity file:- %s\n", RugosityFile);

  strcpy(BVFrequencyFile, tugo_cfg->dissipation->BruntVassalaFile.face_value().c_str()) ;
  fprintf(echo, " Optional BV frequency file:- %s\n", BVFrequencyFile);

  BVFrequencyDefault=tugo_cfg->dissipation->BruntVassalaValue.numerical_value<float>();

/*----------------------------------------------------------------------
  run configuration section */

  strcpy(HotStartFile, tugo_cfg->model->RestartFile.face_value().c_str()) ;

  if(strncmp(HotStartFile, "COLD-START", 10) == 0) {
    fprintf(echo, " Cold start, no file \n");
    modelStart = COLD_START;
    external_state = 0;
    }
  else if(strncmp(HotStartFile, "EXTERNAL", 8) == 0) {
    fprintf(echo, "Hot start from external simulation\n");
    modelStart = EXTERNAL;
    external_state = 1;
    }
  else {
    fprintf(echo, " Hot start file:- %s\n", HotStartFile);
    modelStart = HOT_START;
    external_state = 0;
    }
 
  strcpy(hotstart_formatname, tugo_cfg->model->RestartFormat.face_value().c_str()) ;
  fprintf(echo, "hotstart file format: %s \n", hotstart_formatname);

  hotstart_format = -1;
  if(strcmp(hotstart_formatname, "CLX") == 0)
    hotstart_format = RESTART_FORMAT_ASCII;
  if(strcmp(hotstart_formatname, "NETCDF") == 0)
    hotstart_format = RESTART_FORMAT_NETCDF;
  if(strcmp(hotstart_formatname, "BINARY") == 0)
    hotstart_format = RESTART_FORMAT_BINARY;

  gMinDepth = tugo_cfg->topography->MinDepth.numerical_value<double>() ;
  fprintf(echo, "Absolute minumum depth:- %f\n", gMinDepth);

  gRelativeMinDepth = tugo_cfg->topography->RelativeMinDepth.numerical_value<double>() ;
  fprintf(echo, "Relative minumum depth:- %f\n", gMinDepth);

  gDryingMinDepth = tugo_cfg->topography->DryingMinDepth.numerical_value<double>() ;
  fprintf(echo, "Numerical minumum depth at drying nodes:- %f\n", gDryingMinDepth);

  deltaT[0] = tugo_cfg->model->TimeStep.numerical_value<double>() ;
  fprintf(echo, "Main time step:- %f\n", deltaT[0]);

  nsubcycles=tugo_cfg->misc->nSubCycles.numerical_value<int>();
  if(nsubcycles>=NLVL) {
    check_error(-1,"nsubcycles exceeds max value (NLVL)",__LINE__,__FILE__,1);
    }

//  deltaT[1] = tugo_cfg->model->SubTimeStep.multiple_value<double>((int) 1) ;
/*WARNING*/
//  if(status!=0) deltaT[1]=deltaT[0]/2.0; disabled
//  fprintf(echo, "Secondary time step:- %f\n", deltaT[1]);
  for(k=0;k<nsubcycles;k++) {
    deltaT[k+1] = tugo_cfg->model->SubTimeStep.multiple_value<double>((int) k) ;
    }
  fprintf(echo, "#sub-cycle levels: %d, steps = ", nsubcycles);
  for(k=0;k<=nsubcycles;k++) {
    fprintf(echo, "(level=%d) %e", k,deltaT[k]);
    }
  fprintf(echo, " (seconds) \n");
//  fflush(echo);

  strcpy(pointer, tugo_cfg->model->OrigineTime.face_value().c_str()) ;

  nitems = sscanf(pointer, "%4d", &t_reference.year);

  nitems = sscanf(&pointer[5], "%2d", &t_reference.month) ;

  nitems = sscanf(&pointer[8], "%2d", &t_reference.day) ;

  nitems = sscanf(&pointer[11], "%2d:%2d:%2d", &hour, &Minutes, &Seconds) ;

  t_reference.second = hour*3600. + Minutes*60.+ Seconds;

  fprintf(echo, " Time reference:- day= %d month= %d year =%d\n",
          t_reference.day, t_reference.month, t_reference.year);

  RunDuration = tugo_cfg->model->RunDuration.numerical_value<double>() * 3600. ;

  Days    = (int) (RunDuration / (3600 * 24));
  Hours   = (int) (RunDuration / 3600 - Days * 24);
  Minutes = (int) (RunDuration / 60 - (Days * 24 + Hours) * 60);
  Seconds = (int) (RunDuration - ((Days * 24 + Hours) * 60 + Minutes) * 60);
  fprintf(echo, " Run Length  %d Days,  %d Hours,  %d Minutes,  %f Seconds\n\n",
          Days, Hours, Minutes, Seconds);

  SpinUp = tugo_cfg->model->SpinUpDuration.numerical_value<double>() ;

  fprintf(echo, " Spinup time in hours:- %f\n\n", SpinUp);
  SpinUp *= 3600.;

/*-----------------------------------------------------------------------------
  global tide flag*/
  tidal_forcing = 0;
  
/*-----------------------------------------------------------------------------
  detailed tide flag*/
  LSA_forcing = 0;
  astronomic_forcing=0;
  
  pressure_forcing = 0;
  wind_forcing = 0;
  mean_wind = 0;
  tidal_bc = 0;
  ib_bc = 0;
  spongious_bc = 0;
  time_interpolation = UNDEFINED;

  fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Forcing section\n");

  astronomic_forcing = (tugo_cfg->tides->AstronomicPotentialFlag.actual_value == KEY_TRUE) ;
  if(astronomic_forcing==1) {
    fprintf(echo," Tides forcing included in computation\n");
    }
  else
    fprintf(echo," Tides forcing NOT included in computation\n\n");

  LSA_forcing = (tugo_cfg->tides->LSAFlag.actual_value == KEY_TRUE) ;
  if(LSA_forcing) {
    fprintf(echo," Loading+self-attraction included in computation\n");
    }
  else
    fprintf(echo," Loading+self-attraction NOT included in computation\n\n");

  SpectralAtmosphericPressure = (tugo_cfg->tides->PressureFlag.actual_value == KEY_TRUE) ;
  if(SpectralAtmosphericPressure) {
    fprintf(echo," Spectral atmospheric surface pressure forcing included in computation\n");
    }
  else
    fprintf(echo," Spectral atmospheric surface pressure NOT included in computation\n\n");

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:

  01/09/2008 :

    a better method is needed to enable/disable the atmospheric
    contributions in computation

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  if(tugo_cfg->atmosphere->OnOffFlag.actual_value == KEY_FALSE) {
    fprintf(echo," Atmospheric effects globally disabled in computation\n");
    tugo_cfg->atmosphere->PressureFlag.actual_value=KEY_FALSE;
    tugo_cfg->atmosphere->WindFlag.actual_value=KEY_FALSE;
    tugo_cfg->atmosphere->KeepMeanWindFlag.actual_value=KEY_FALSE;
    tugo_cfg->atmosphere->BoundInversebarometer.actual_value=KEY_FALSE;
    }
  else
    fprintf(echo, " Atmospheric effects enabled in computation\n");

  pressure_forcing = (tugo_cfg->atmosphere->PressureFlag.actual_value == KEY_TRUE) ;

  if(pressure_forcing) {
    fprintf(echo, " Pressure forcing included in computation\n");
    }
  else
    fprintf(echo, " Pressure forcing NOT included in computation\n");

  wind_forcing = (tugo_cfg->atmosphere->WindFlag.actual_value == KEY_TRUE) ;

  if(wind_forcing) {
    fprintf(echo, " Wind forcing included in computation\n");
    }
  else
    fprintf(echo, " Wind forcing NOT included in computation\n\n");

  mean_wind = (tugo_cfg->atmosphere->KeepMeanWindFlag.actual_value == KEY_TRUE) ;

  if(mean_wind) {
    fprintf(echo, " Do not remove mean wind from forcing\n");
    }
  else
    fprintf(echo, "Remove mean wind from forcing\n\n");

  time_interpolation = UNDEFINED;
  strcpy(key, tugo_cfg->atmosphere->TimeInterp.face_value().c_str()) ;
  if (strcmp(key, "LINEAR") == 0) {
    fprintf(echo, " Forcing interpolation is LINEAR in time\n");
    time_interpolation = LINEAR;
    }
    
  if (strcmp(key, "QUADRATIC") == 0) {
    fprintf(echo, " Forcing interpolation is QUADRATIC in time\n\n");
    time_interpolation = QUADRATIC;
    };
/*
  if(strcmp(key,"CUBIC")==0) {
    fprintf(echo," Forcing interpolation is CUBIC in time\n\n");
    time_interpolation=CUBIC;
    };
*/
  if(time_interpolation == UNDEFINED) {
    check_error(-1, "Forcing interpolation input is incorrect", __LINE__,__FILE__, 1);
    };

  bottom_forcing = tugo_cfg->topography->BottomMotionFlag.numerical_value() ;
  if(bottom_forcing) {
    fprintf(echo, " Bottom motion forcing included in computation\n");
    }
  else
    fprintf(echo, " Bottom motion forcing NOT included in computation\n\n");

  ocean_loading = tugo_cfg->misc->NonTidalLoadFlag.numerical_value() ;
  if(ocean_loading) {
    fprintf(echo, " Non tidal ocean loading forcing included in computation\n");
    }
  else
    fprintf(echo, "  Non tidal ocean forcing NOT included in computation\n\n");

  g_oceanwaves_forcing = (tugo_cfg->OceanWaves->OnOffFlag.actual_value == KEY_TRUE) ;
  g_oceanwaves_time_interpolation = UNDEFINED;
  if(g_oceanwaves_forcing) {
    fprintf(echo, " OceanWaves included in computation\n");
    strcpy(key, tugo_cfg->OceanWaves->TimeInterp.face_value().c_str()) ;
    if (strcmp(key, "LINEAR") == 0) {
      fprintf(echo, " OceanWaves Forcing interpolation is LINEAR in time\n");
      g_oceanwaves_time_interpolation = LINEAR;
      };
    if (strcmp(key, "QUADRATIC") == 0) {
      fprintf(echo, " Oceanwaves Forcing interpolation is QUADRATIC in time\n\n");
      g_oceanwaves_time_interpolation = QUADRATIC;
      };
    if(g_oceanwaves_time_interpolation == UNDEFINED) {
      check_error(-1, "OceanWaves Forcing interpolation input is incorrect", __LINE__,__FILE__, 1);
      };
    }
  else
    fprintf(echo, " OceanWaves NOT included in computation\n\n");

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

/*----------------------------------------------------------------------
  Boundary conditions section */

  fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx BCs section\n");
  ElevationOBC_type = tugo_cfg->boundaries->type.numerical_value() ;
  fprintf(echo, " BC's routine : %d\n", ElevationOBC_type);

  if(tugo_cfg->boundaries->mode.face_value() == "LOCAL") {
    ElevationOBC_mode=BOUNDARY_LOCAL;
    }
  else if (tugo_cfg->boundaries->mode.face_value() == "VARIATIONAL") {
    ElevationOBC_mode=BOUNDARY_VARIATIONAL;
    }
  else {
    sprintf(msg, "unreckognized option for boundaries->mode : %s\n",tugo_cfg->boundaries->mode.face_value().c_str());
    check_error(-1, msg, __LINE__, __FILE__, 1);
    }

  fprintf(echo, " BC's routine : %d\n", ElevationOBC_type);

  tidal_bc = (tugo_cfg->tides->BoundTideFlag.actual_value == KEY_TRUE) ;
  if(tidal_bc) {
    fprintf(echo, " Tides BC's included in computation\n");
    }
  else
    fprintf(echo, " Tides BC's NOT included in computation\n");

  ib_bc = (tugo_cfg->atmosphere->BoundInversebarometer.actual_value == KEY_TRUE) ;
  if(ib_bc) {
    fprintf(echo, " Inverse barometer BC's included in computation\n");
    }
  else
    fprintf(echo, " Inverse barometer BC's NOT included in computation\n");

  external_obc = (tugo_cfg->boundaries->ArchivedBCUse.actual_value == KEY_TRUE) ;
  if(external_obc) {
    fprintf(echo, " External archive BC's included in computation\n");
    }
  else
    fprintf(echo, " External BC's NOT included in computation\n");

  spongious_bc = (tugo_cfg->boundaries->SpongeFlag.actual_value == KEY_TRUE) ;
  if(spongious_bc) {
    fprintf(echo, " Spongious buffer at the open limits\n");
    }
  else
    fprintf(echo, " NO spongious buffer at the open limits\n");

/*----------------------------------------------------------------------
  ??? section */

  smooth_instabilities = (tugo_cfg->misc->SmoothingFlag.actual_value == KEY_TRUE) ;
  if(smooth_instabilities) {
    fprintf(echo, " Sea state vector smoothed when instable\n");
    }
  else
    fprintf(echo, " Sea state vector NOT smoothed when instable\n");

  SmoothZField = (tugo_cfg->misc->SmoothZField.actual_value == KEY_TRUE) ;
  if(SmoothZField) {
    fprintf(echo, " Z field smoothed (sort of diffusion)\n");
    }
  else
    fprintf(echo, " Z field NOT smoothed\n");

  SmoothUField = (tugo_cfg->misc->SmoothUField.actual_value == KEY_TRUE) ;
  if(SmoothUField) {
    fprintf(echo, " U field smoothed (sort of diffusion)\n");
    }
  else
    fprintf(echo, " U field NOT smoothed\n");

  SmoothAdvection = (tugo_cfg->misc->SmoothAdvection.actual_value == KEY_TRUE) ;
  if(SmoothAdvection) {
    fprintf(echo, " Advection smoothed (sort of diffusion)\n");
    }
  else
    fprintf(echo, " Advection NOT smoothed\n");

  SmoothDiffusion = (tugo_cfg->misc->SmoothDiffusion.actual_value == KEY_TRUE) ;
  if(SmoothDiffusion) {
    fprintf(echo, " Diffusion smoothed (sort of diffusion)\n");
    }
  else
    fprintf(echo, " Diffusion NOT smoothed\n");

  Spectral2DRun = (tugo_cfg->tides->Spectral2DRun.actual_value == KEY_TRUE) ;
  if(Spectral2DRun) {
    fprintf(echo, " Solve 2D-spectral tides\n");
    }
  else
    fprintf(echo, " Do NOT solve 2D-spectral tides\n");

  if(     tugo_cfg->tides->Spectral3DRun.actual_value == (string) "BAROTROPIC") {
    fprintf(echo, " Solve 3D-spectral barotropic tides\n");
    Spectral3DRun =BAROTROPIC_MODE;
    }
  else if(tugo_cfg->tides->Spectral3DRun.actual_value == (string) "BAROCLINIC") {
    fprintf(echo, " Solve 3D-spectral baroclinic tides\n");
    Spectral3DRun =BAROCLINIC_MODE;
    }
  else if(tugo_cfg->tides->Spectral3DRun.actual_value == (string) "SW2D+1D") {
    fprintf(echo, " Solve 2D-spectral + 1D barotropic tides\n");
    Spectral3DRun =SW2D1D_MODE;
    }
  else {
    fprintf(echo, " Do NOT solve 3D-spectral tides\n");
    Spectral3DRun =UNDEFINED_MODE;
    }

  admittance = (tugo_cfg->tides->admittance.actual_value == KEY_TRUE) ;
  if(admittance) {
    fprintf(echo, " Admittance extension for tidal spectrum \n");
    }
  else
    fprintf(echo, " NO admittance extension for tidalspectrum\n");

  nodal_corrections = (tugo_cfg->tides->NodalCorrection.actual_value == KEY_TRUE) ;
  if(nodal_corrections) {
    fprintf(echo, " Nodal correction included in computation and analysis\n");
    }
  else
    fprintf(echo, " Nodal correction NOT included in computation NOR in analysis\n");

/*----------------------------------------------------------------------
  global tidal atlas directory*/
  strcpy(atlas_directory, tugo_cfg->tides->AtlasDir.face_value().c_str()) ;
  fprintf(echo, "tidal atlas directory name: %s\n", atlas_directory);

/*----------------------------------------------------------------------
  global tidal atlas finename convention*/
  strcpy(atlas_convention, tugo_cfg->tides->AtlasConvention.face_value().c_str()) ;
  fprintf(echo, "tidal atlas name convention: %s\n", atlas_convention);

  need_meteo = (pressure_forcing || wind_forcing || ib_bc);

  fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Solver section\n");
  strcpy(solver_name, tugo_cfg->model->LinearSolver.face_value().c_str()) ;
  fprintf(echo, "main solver name : %s \n", solver_name);

  gLinearSolverID =  LinearSystem_identify(solver_name,false);
//   gLinearSolverID =  LinearSystem_identify(solver_name);
  printf("solver name %s internal reference = %d\n", solver_name, gLinearSolverID);
  if(gLinearSolverID==-1) {
    check_error(-1, "prescribed solver not implemented", __LINE__, __FILE__, 1);
    }
    
  strcpy(tmp, tugo_cfg->model->SolverMode.face_value().c_str()) ;
  fprintf(echo, " solver mode : %s \n", tmp);
  
  gLinearSolverMode=-1;
  
  if(strcmp(tmp, "AUTO")   == 0)
    gLinearSolverMode=-1;
  if(strcmp(tmp, "SEQUENTIAL") == 0)
    gLinearSolverMode=SEQUENTIAL_COMPUTING;
  if(strcmp(tmp, "PARALLEL") == 0)
    gLinearSolverMode=PARALLEL_COMPUTING;
  

  strcpy(solver_name, tugo_cfg->model->SubLinearSolver.face_value().c_str()) ;
  fprintf(echo, "subcycling solver name : %s \n", solver_name);

  if(strcmp(solver_name, "AUTO") == 0)
    gSubLinearSolverID = gLinearSolverID;
  else
    gSubLinearSolverID =  LinearSystem_identify(solver_name,false);
//     gSubLinearSolverID =  LinearSystem_identify(solver_name);
  printf("subcycling solver name %s internal reference = %d\n", solver_name, gSubLinearSolverID);
  if(gSubLinearSolverID==-1) {
    check_error(-1, "prescribed solver not implemented", __LINE__, __FILE__, 1);
    }

  strcpy(solver_name, tugo_cfg->tides->SpectralSolver.face_value().c_str()) ;
  fprintf(echo, "spectral solver name : %s \n", solver_name);

  if(strcmp(solver_name, "AUTO") == 0)
    gSpectralSolver = gLinearSolverID;
  else
    gSpectralSolver =  LinearSystem_identify(solver_name,false);
//     gSpectralSolver =  LinearSystem_identify(solver_name);
  printf("spectral solver name %s internal reference = %d\n", solver_name, gSpectralSolver);
  if(gSpectralSolver==-1) {
    check_error(-1, "prescribed slover not implemented", __LINE__, __FILE__, 1);
    }

/**----------------------------------------------------------------------------
  Shallow-water solver */
  setting=tugo_cfg->dynamics2D->WaveEquationForm.face_value();

  if (setting      == KEY_GWE_CLX) {
    we_type = GWE_CLX ;
    gFormulation=FINITE_ELEMENTS;
    }
  else if (setting == KEY_GWE_XXX) {
    we_type = GWE_XXX ;
    gFormulation=FINITE_ELEMENTS;
    }
  else if (setting == KEY_WE_EXPLICIT) {
    we_type = WE_EXPLICIT ;
    gFormulation=FINITE_ELEMENTS;
    }
  else if (setting == KEY_WE_IMPLICIT) {
    we_type = WE_IMPLICIT ;
    gFormulation=FINITE_ELEMENTS;
    }
  else if (setting == KEY_FV_IMPLICIT) {
    we_type = FV_IMPLICIT ;
    gFormulation=FINITE_VOLUMES;
    }
  else {
    printf("unknown shallow-water solver: %s\n",setting.c_str());
    check_error(-1, "unknown shallow-water solver", __LINE__, __FILE__, 1);
    }

  fprintf(echo, "Wave equation mode:- %d\n", we_type);

/*----------------------------------------------------------------------
  physical parameters section */

  explicit_loading = tugo_cfg->tides->LSACoeff.numerical_value<double>() ;
  fprintf(echo, " Explicit tidal LSA coefficient:- %f\n", explicit_loading);

  if((explicit_loading != 1.0) && (ocean_loading == 0)) {
    fprintf(echo,"WARNING: ocean_loading not used and explicit tidal LSA coefficient not equal to 1\n\n");
    }

  linear_h = tugo_cfg->tides->DeformationLoveNum.numerical_value<double>() ;
  fprintf(echo, " Linear deformation Love number:- %lf\n", linear_h);

  fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Dissipation section\n");
  
  gR_default = tugo_cfg->dissipation->LinearFrictionCoeff.numerical_value<double>() ;
  fprintf(echo, " Default linear drag coefficient:- %f\n", gR_default);

  gR3D_default = tugo_cfg->dissipation3D->LinearFrictionCoeff.numerical_value<double>() ;
  fprintf(echo, " Default 3D linear drag coefficient:- %f\n", gR3D_default);

  gCd_default = tugo_cfg->dissipation->QuadraticFrictionCoeff.numerical_value<double>() ;
  fprintf(echo, " Default quadratic drag coefficient:- %f\n", gCd_default);

  gCd_minimum=tugo_cfg->dissipation->MinQuadraticFrictionCoeff.numerical_value<double>();
  fprintf(echo, " Default minimum quadratic drag coefficient:- %f\n", gCd_minimum);
  
  gz0_default=tugo_cfg->dissipation->RoughnessLength.numerical_value<double>();
  fprintf(echo, " Default 2D rugosity length :- %f\n", gz0_default);
  
  gz03D_default=tugo_cfg->dissipation3D->RoughnessLength.numerical_value<double>();
  fprintf(echo, " Default 3D rugosity length :- %f\n", gz03D_default);
  
  if (tugo_cfg->dissipation->BottomFrictionType.face_value() == KEY_LINEAR) {
    gBottomFrictionType=FRICTION_LINEAR;
    fprintf(echo, " Bottom friction parameterisation:- %s\n", KEY_LINEAR.c_str());
    }
  else if (tugo_cfg->dissipation->BottomFrictionType.face_value() == KEY_QUADRATIC) {
/**----------------------------------------------------------------------
    Quadratic friction with Cd formally prescribed*/
    gBottomFrictionType=FRICTION_QUADRATIC;
    fprintf(echo, " Bottom friction parameterisation:- %s\n", KEY_QUADRATIC.c_str());
    }
  else if (tugo_cfg->dissipation->BottomFrictionType.face_value() == KEY_KARMAN) {
/**----------------------------------------------------------------------
    Quadratic friction with Cd deduced from velocity log-profile assumption*/
    gBottomFrictionType=FRICTION_KARMAN;
    fprintf(echo, " Bottom friction parameterisation:- %s\n", KEY_KARMAN.c_str());
    }
  else {
    check_error(-1, "Invalid bottom friction type", __LINE__, __FILE__, 1) ;
    }

  gCm = tugo_cfg->dissipation->MixedLayerCoeff.numerical_value<double>() ;
  fprintf(echo, " Mixed layer/shear drag coefficient :- %f\n", gCm);

  gU0 = tugo_cfg->dissipation->MinBackGroundSpd.numerical_value<double>() ;
  if(gU0>1.0) {
    printf("!!! Minimum background friction speed %f apparently given in cm/s (i.e. u0 > 1m/s) , pleas check\n", gU0);
    fprintf(echo, "!!! Minimum background friction speed %f apparently given in cm/s (i.e. u0 > 1m/s) , pleas check\n", gU0);
//     fprintf(echo,"!!! Minimum background friction speed %f apparently given in cm/s, changed to m/s\n", gU0);
//     gU0/=100.;
    }
  fprintf(echo, " Minimum background friction speed (m/s):- %f\n", gU0);

  minKh = tugo_cfg->dissipation->MinhorizontalVisc.numerical_value<double>() ;
  if(minKh>100.0) {
    printf("!!! Minimum horizontal eddy viscosity %f apparently given in cm²/s, changed to m²/s\n", minKh);
    fprintf(echo,"!!! Minimum horizontal eddy viscosity %f apparently given in cm²/s, changed to m²/s\n", minKh);
    minKh/=10000.;
    }
  fprintf(echo, " Minimum horizontal eddy viscosity (m^2/s):- %f\n", minKh);

  tau0 = tugo_cfg->dynamics2D->tau0.numerical_value<double>() ;
  fprintf(echo, " Tau0:- %f\n", tau0);

  gTheta = tugo_cfg->dynamics2D->theta.numerical_value<double>() ;
#ifdef DRY
  /* Fred's dry change: explicit timestepping */
  gTheta = 0.0;
#endif
  fprintf(echo, " Theta:- %f\n\n", gTheta);

  tau_r = tugo_cfg->boundaries->RelaxTime.numerical_value<double>() ;
  fprintf(echo, " Relaxation (hours):- %f\n\n", tau_r);
  tau_r *= 3600.;

/*----------------------------------------------------------------------
  runtime IO section */
  fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Runtime IO section\n");
  
  strcpy(TideDataFile, tugo_cfg->tides->BoundTideFile.face_value().c_str()) ;
  fprintf(echo, "tidal boundary conditions file: %s\n", TideDataFile);

  strcpy(external_path, tugo_cfg->boundaries->ArchivedBCDirectory.face_value().c_str()) ;
  fprintf(echo, "external OBC directory name: %s \n", external_path);

  strcpy(external_root, tugo_cfg->boundaries->ArchivedBCRoot.face_value().c_str()) ;
  fprintf(echo, "external OBC file root name: %s \n", external_root);

  strcpy(LSA_directory, tugo_cfg->tides->LSADir.face_value().c_str()) ;
  fprintf(echo, "Loading/self-attraction potential directory name: %s\n", LSA_directory);

  strcpy(LSA_convention, tugo_cfg->tides->LSAConvention.face_value().c_str()) ;
  fprintf(echo, "Loading/self-attraction name convention: %s\n", LSA_convention);

  strcpy(SAP_directory, tugo_cfg->tides->PressureDir.face_value().c_str()) ;
  fprintf(echo, "spectral atmospheric surface pressure directory name: %s\n", SAP_directory);

  strcpy(SAP_convention, tugo_cfg->tides->PressureConvention.face_value().c_str()) ;
  fprintf(echo, "spectral atmospheric surface pressure convention: %s\n", SAP_convention);

  strcpy(ForcingD, tugo_cfg->atmosphere->AtmosDirectory.face_value().c_str()) ;
  fprintf(echo, "atmospheric forcing directory name: %s \n", ForcingD);

  strcpy(meteo_formatname, tugo_cfg->atmosphere->FileFormat.face_value().c_str()) ;
  fprintf(echo, "atmospheric forcing file format: %s \n", meteo_formatname);

  meteo_format = -1;
  if(strcmp(meteo_formatname, "CLX") == 0)
    meteo_format = 0;
  if(strcmp(meteo_formatname, "NETCDF") == 0)
    meteo_format = 1;
  if(strcmp(meteo_formatname, "ACADEMIC") == 0)
    meteo_format = 2;
  if(strcmp(meteo_formatname, "HARMONIC") == 0)
    meteo_format = 3;
  if(strcmp(meteo_formatname, "GRIB") == 0)
    meteo_format = 4;

  strcpy(meteo_convention, tugo_cfg->atmosphere->FileConvention.face_value().c_str()) ;
  fprintf(echo, "atmospheric forcing naming convention: %s \n", meteo_convention);

  strcpy(meteo_maskfile, tugo_cfg->atmosphere->MaskFile.face_value().c_str()) ;
  fprintf(echo, "atmospheric forcing mask: %s \n",meteo_maskfile);
  
/*---------------------------------------------------------------------------------
/* Runtime I/O OceanWaves section **/
  
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
* \file configure-hugo.cpp
* \date 26/01/2010
*       \todo Clément : could we do the same by reading directly the wave format in tugo_cfg when we need it
*              instead of using th global variable waves_format ?
*
**/

  if  (tugo_cfg->OceanWaves->OnOffFlag.actual_value == KEY_TRUE) {
    strcpy(gWavesDirectory, tugo_cfg->OceanWaves->FileDirectory.face_value().c_str()) ;
    fprintf(echo, "atmospheric forcing directory name: %s \n", ForcingD);
    strcpy(g_waves_formatname, tugo_cfg->OceanWaves->FileFormat.face_value().c_str()) ;
    fprintf(echo, "ocean waves forcing file format: %s \n", g_waves_formatname);

    g_waves_format = -1;
    if(strcmp(g_waves_formatname, "WW3_SG_NETCDF") == 0)
      g_waves_format = WW3_SG_NETCDF;                       //pre-processor : define WW3_SG_NETCDF 0
    if(strcmp(g_waves_formatname, "WW3_UG_NETCDF") == 0)
      g_waves_format = WW3_UG_NETCDF;                       //pre-processor : define WW3_UG_NETCDF 1
    if(strcmp(g_waves_formatname, "WW3_SG_NATIVE") == 0)
      g_waves_format = WW3_SG_NATIVE;                       //pre-processor : define WW3_SG_NETCDF 2
    if(strcmp(g_waves_formatname, "WW3_UG_NATIVE") == 0)
      g_waves_format = WW3_UG_NATIVE;                       //pre-processor : define WW3_UG_NETCDF 3

    strcpy(gWavesConvention, tugo_cfg->OceanWaves->FileConvention.face_value().c_str()) ;
    fprintf(echo, "ocean waves forcing naming convention: %s \n", gWavesConvention);

    strcpy(g_waves_maskfile, tugo_cfg->OceanWaves->MaskFile.face_value().c_str()) ;
    fprintf(echo, "atmospheric forcing mask: %s \n",g_waves_maskfile);
    }
  
  /*----------------------------------------------------------------------------------*/

  strcpy(quaker_directory, tugo_cfg->topography->BottomMotionDirectory.face_value().c_str()) ;
  fprintf(echo, "bottom motion directory name: %s \n", quaker_directory);

  strcpy(quaker_convention, tugo_cfg->topography->BottomMotionFileConvent.face_value().c_str()) ;
  fprintf(echo, "bottom motion naming convention: %s \n", quaker_convention);

  strcpy(ContinueFile, tugo_cfg->model->ContinuationFile.face_value().c_str()) ;
  fprintf(echo, "rootname for continuation files: %s \n", ContinueFile);

  strcpy(restart_formatname, tugo_cfg->model->ContinuationFormat.face_value().c_str()) ;
  fprintf(echo, "restart file format: %s \n", restart_formatname);

  restart_format = -1;
  if(strcmp(restart_formatname, "CLX") == 0)
    restart_format = RESTART_FORMAT_ASCII;
  if(strcmp(restart_formatname, "NETCDF") == 0)
    restart_format = RESTART_FORMAT_NETCDF;
  if(strcmp(restart_formatname, "BINARY") == 0)
    restart_format = RESTART_FORMAT_NETCDF;

  ContFInt = tugo_cfg->model->ContinuationInterval.numerical_value<double>() ;
  Days    = (int) (ContFInt / (3600 * 24));
  Hours   = (int) (ContFInt / 3600 - Days * 24);
  Minutes = (int) (ContFInt / 60 - (Days * 24 + Hours) * 60);
  Seconds = (int) (ContFInt - ((Days * 24 + Hours) * 60 + Minutes) * 60);
  fprintf(echo, "Written every %d Days, %d Hours,  %d Minutes, %f Seconds\n",
          Days, Hours, Minutes, Seconds);

/*----------------------------------------------------------------------
  unstructured grid archives */
  fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx FE archive section\n");
  if (tugo_cfg->fearchive->OnOffFlag.actual_value == KEY_TRUE) {
    fprintf(echo, " Runtime saving ON\n");
    runtime_save=1;
    }
  else {
    fprintf(echo, " Runtime saving OFF\n");
    runtime_save=0;
    }

  strcpy(UGarchive_rootname, tugo_cfg->fearchive->StateFile.face_value().c_str()) ;
  fprintf(echo, "name for state vector binary file: %s \n", UGarchive_rootname);

  strcpy(UGforcing_rootname, tugo_cfg->fearchive->VectorFile.face_value().c_str()) ;
  fprintf(echo, "name for forcing vector binary file: %s \n", UGforcing_rootname);

  strcpy(archive_formatname, tugo_cfg->fearchive->ArchFormat.face_value().c_str()) ;
  fprintf(echo, "archive file format: %s \n", archive_formatname);

  archive_format = -1;
  if(strcmp(archive_formatname, "CLX") == 0)
    archive_format = 0;
  if(strcmp(archive_formatname, "NETCDF") == 0)
    archive_format = 1;
  if(strcmp(archive_formatname, "NETCDF-OBSOLETE") == 0)
    archive_format = 2;

  strcpy(archive_convention, tugo_cfg->fearchive->ArchConvention.face_value().c_str()) ;
  fprintf(echo, "archive naming convention: %s \n", archive_convention);

  archive_start = tugo_cfg->fearchive->ArchStart.numerical_value<double>() ;
  fprintf(echo, "State vector and forcing storage start   : %lf (days)\n",archive_start);
  archive_start = archive_start * 24. * 3600.;

  archive_sampling = tugo_cfg->fearchive->ArchInterval.numerical_value<double>() ;
  fprintf(echo, "State vector and forcing storage interval: %lf (seconds)\n", archive_sampling);

  archive_type = tugo_cfg->fearchive->ArchBits.numerical_value() ;
  fprintf(echo, "State vector and forcing storage type: %d (bytes)\n",archive_type);

  sscanf(tugo_cfg->fearchive->Arch_huv_Scale.actual_value.c_str(),
         "%lf %lf %lf", &archive_scale1[0], &archive_scale1[1], &archive_scale1[2]) ;
  sscanf(tugo_cfg->fearchive->Arch_pW_Scale.actual_value.c_str(),
         "%lf %lf %lf", &archive_scale2[0], &archive_scale2[1],&archive_scale2[2]) ;
  fprintf(echo, " Archive scale 1: %f %f %f \n", archive_scale1[0],archive_scale1[1], archive_scale1[2]);
  fprintf(echo, " Archive scale 2: %f %f %f \n", archive_scale2[0],archive_scale2[1], archive_scale2[2]);

  if( tugo_cfg->fearchive->OverwriteAppendFlag.face_value() == "OVERWRITE") {
    archive_update = 0;
    }
  else {
    archive_update = 1;
    }
  fprintf(echo, " Allow for existing archives update %d (flag)\n",archive_update);

/*----------------------------------------------------------------------
  structured grid archives */
  fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx SG archive section\n");
  if (tugo_cfg->structured->OnOffFlag.actual_value == Keywords::KEY_TRUE) {
    regular_output=1;
    fprintf(echo, " Structured grid archive ON\n");
    }
  else {
    regular_output=0;
    fprintf(echo, " Structured grid archive OFF\n");
    }

  strcpy(gStructuredZone, tugo_cfg->structured->zone.face_value().c_str()) ;
  fprintf(echo, "name for structured grid archive file: %s \n", gStructuredZone);

  regular_start = tugo_cfg->structured->Start.numerical_value<double>() ;
  fprintf(echo, "Structured grid archive storage start   : %lf (days)\n",regular_start);
  regular_start = regular_start * 24. * 3600.;

  regular_sampling = tugo_cfg->structured->Interval.numerical_value<double>() ;
  fprintf(echo, "Structured grid archive storage interval: %lf (seconds)\n", regular_sampling);

/*----------------------------------------------------------------------
  local sampling archives */
  if (tugo_cfg->sample->OnOffFlag.actual_value == KEY_TRUE) {
    sampling=1;
    fprintf(echo, " Sample saving ON\n");
    }
  else {
    sampling=0;
    fprintf(echo, " Sample saving OFF\n");
    }

  strcpy(PlotInputFile, tugo_cfg->sample->file.face_value().c_str()) ;
  fprintf(echo, "name for sample points input file: %s \n", PlotInputFile);

  PlotStart = tugo_cfg->sample->start.numerical_value<double>() ;
  fprintf(echo, "Sample storage start   : %lf \n", PlotStart);

  PlotInterval = tugo_cfg->sample->SaveInterval.numerical_value<double>() ;
  fprintf(echo, "Sample storage interval: %lf \n", PlotInterval);

  if( tugo_cfg->sample->OverAppendFlag.face_value()== (string) "OVERWRITE") {
    plot_update = 0;
    }
  else if( tugo_cfg->sample->OverAppendFlag.face_value()== (string) "APPEND") {
    plot_update = 1;
    }
//   if( strcmp(tugo_cfg->sample->OverAppendFlag.face_value().c_str(),"OVERWRITE")==0) {
//     plot_update = 0;
//     }
//   else if( strcmp(tugo_cfg->sample->OverAppendFlag.face_value().c_str(),"APPEND")==0) {
//     plot_update = 1;
//     }
  else {
    sprintf(msg, "unreckognized option for sample->OverAppendFlag (OVERWRITE/APPEND) : %s\n",tugo_cfg->sample->OverAppendFlag.face_value().c_str());
    check_error(-1, msg, __LINE__, __FILE__, 1);
    }
//  plot_update = tugo_cfg->sample->OverAppendFlag.numerical_value() ;
  fprintf(echo, " Allow for existing sampling update %d (flag)\n", plot_update);

  strcpy(SectionInputFile, tugo_cfg->sample->SectionInputFile.face_value().c_str()) ;
  fprintf(echo, "name for sample sections input file: %s \n", SectionInputFile);

/*----------------------------------------------------------------------
  harmonic analysis section */
  fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Harmonic analysis section\n");
  if (tugo_cfg->analysis->OnOffFlag.actual_value == KEY_TRUE) {
    harmonic=1;
    fprintf(echo, " Harmonic analysis ON\n");
    }
  else {
    harmonic=0;
    fprintf(echo, " Harmonic analysis OFF\n");
    }

  gHarmonicAcquisitionStart = tugo_cfg->analysis->HarmoStart.numerical_value<double>() ;
  fprintf(echo, "Harmonic storage start : %lf (days)\n", gHarmonicAcquisitionStart);
  gHarmonicAcquisitionStart = gHarmonicAcquisitionStart * 24. * 3600.;

  gHarmonicAcquisitionInterval = tugo_cfg->analysis->SampleInterval.numerical_value<double>() ;
  fprintf(echo, "Harmonic storage interval : %lf\n", gHarmonicAcquisitionInterval);

  gHarmonicAnalysisInterval = tugo_cfg->analysis->ComputeInterval.numerical_value<double>() ;
  fprintf(echo, "Harmonic analysis interval: %lf (days)\n", gHarmonicAnalysisInterval);
  gHarmonicAnalysisInterval = gHarmonicAnalysisInterval * 24. * 3600.;

  if (tugo_cfg->analysis->SequentialRecycle.actual_value == (string) "MANAGED") {
    gHarmonic_RecycleSequential=-1;
    fprintf(echo, " Harmonic solution recycling for OBCs MANAGED\n");
    }
  else if (tugo_cfg->analysis->SequentialRecycle.actual_value == KEY_TRUE) {
    gHarmonic_RecycleSequential=1;
    fprintf(echo, " Harmonic solution recycling for OBCs ON\n");
    }
  else {
    gHarmonic_RecycleSequential=0;
    fprintf(echo, " Harmonic solution recycling for OBCs OFF\n");
    }

  if (tugo_cfg->analysis->SGArchiveFlag.actual_value == KEY_TRUE) {
    gHarmonic_StructuredArchive=1;
    fprintf(echo, " Harmonic solution archived on regular grid ON\n");
    }
  else {
    gHarmonic_StructuredArchive=0;
    fprintf(echo, " Harmonic solution archived on regular grid OFF\n");
    }

  if (tugo_cfg->analysis->LGP1ArchiveFlag.actual_value == KEY_TRUE) {
    gHarmonic_LGP1Archive=1;
    fprintf(echo, " Harmonic solution archived on LGP1 discretisation (ASCII) ON\n");
    }
  else {
    gHarmonic_LGP1Archive=0;
    fprintf(echo, " Harmonic solution archived on LGP1 discretisation (ASCII)  OFF\n");
    }

  if (tugo_cfg->analysis->NETCDFArchiveFlag.actual_value == KEY_TRUE) {
    gHarmonic_NETCDFArchive=1;
    fprintf(echo, " Harmonic solution archived on NETCDF format ON\n");
    }
  else {
    gHarmonic_NETCDFArchive=0;
    fprintf(echo, " Harmonic solution archived on NETCDF format  OFF\n");
    }

  need_tides = (astronomic_forcing || LSA_forcing || tidal_bc || harmonic );

/*----------------------------------------------------------------------
  energy budget section */
  energetics = 0;
  gSpectralEnergy  = 0;

  if(tugo_cfg->energy->OnOffFlag.actual_value == KEY_TRUE) {
    fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Energy processing section\n");
    energetics = (tugo_cfg->energy->TimeStepComp.actual_value == KEY_TRUE) ;
    if(energetics) {
      fprintf(echo, " Energy checks ON\n");
      TEnergyFile=gOutputPath + (string) "/TEnergyFile";
      }
    else
      fprintf(echo, " Energy checks OFF\n");

    EnergyStart = tugo_cfg->energy->start.numerical_value<double>() ;
    fprintf(echo, "Energy storage start   : %lf (seconds)\n", EnergyStart);

    EnergyInterval = tugo_cfg->energy->interval.numerical_value<double>() ;
    fprintf(echo, "Energy storage interval: %lf (seconds)\n", EnergyInterval);

/**----------------------------------------------------------------------------
    MPI unsafe */
    if(energetics) {
      out_budget = fopen("budget.out", "w");
      }

/*----------------------------------------------------------------------
    spectral energy budget section */
    if (tugo_cfg->energy->HarmoComp.actual_value == KEY_TRUE) {
      gSpectralEnergy=1;
      fprintf(echo, " spectral energy budget ON\n");
      }
    else {
      gSpectralEnergy=0;
      fprintf(echo, " spectral energy budget OFF\n");
      }
    }


/*----------------------------------------------------------------------
  additional check section */
  checks = tugo_cfg->misc->checks.numerical_value() ;
  if(checks) {
    fprintf(echo, "Miscallenous checks ON\n");
    }
  else
    fprintf(echo, " Miscallenous checks OFF\n");

/**----------------------------------------------------------------------------
  MPI unsafe */
  if(checks) {
    out_checks = fopen("checks.out", "w");
    }

/*----------------------------------------------------------------------
  UGO section */
  if (tugo_cfg->model->ThreeDPressure.actual_value == KEY_TRUE) {
    mean_3Dpressure=1;
    fprintf(echo, " Mean baroclinic pressure included in computation\n");
    }
  else {
    mean_3Dpressure=0;
    fprintf(echo, " Mean baroclinic pressure NOT included in computation\n\n");
    }

  compatibility_scale=1.0;
  compatibility_scale = tugo_cfg->misc->CompatibilityScale.numerical_value<double>() ;
  fprintf(echo, "Compatibility scale   : %lf\n", compatibility_scale);

  /*----------------------------------------------------------------------
  ICE section */

  if (tugo_cfg->ice->OnOffFlag.actual_value == KEY_TRUE) {
    use_ice=1;
    fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Sea ice section\n");
    fprintf(echo," Ice module ON\n");

    sscanf(tugo_cfg->ice->IceCoverFile.face_value().c_str(), "%s ", icefile) ;

    fprintf(echo,"icefile for concentration: %s\n",icefile);
    dtice = tugo_cfg->ice->TimeStep.numerical_value<double>() ;
    fprintf(echo,"timestep for ice module  : %f (seconds)\n",dtice);
    ice_mom = tugo_cfg->ice->MomentumSolver.numerical_value() ;
    switch (ice_mom) {
      case (0):
/*----------------------------------------------------------------------
        free-drift */
        fprintf(echo,"ice momemtum solver      : free-drift\n");
        break;

      case (1):
/*----------------------------------------------------------------------
        Viscous-Plastic */
        fprintf(echo,"ice momemtum solver      : viscous-plastic\n");
        break;

      case (2):
/*----------------------------------------------------------------------
        EVP */
        fprintf(echo,"ice momemtum solver      : EVP\n");
        break;

      case (3):
/*----------------------------------------------------------------------
        Cavitating */
        fprintf(echo,"ice momemtum solver      : cavitating\n");
        break;
      }

    ice_rheo= tugo_cfg->ice->Rheology.numerical_value() ;
    switch (ice_rheo) {
      case (0):
/*----------------------------------------------------------------------
        VP rheology */
        fprintf(echo,"ice rheology             : VP\n");
        break;

      case (1):
/*----------------------------------------------------------------------
        Mohr-Coulomb */
        fprintf(echo,"ice rheology             : Mohr-Coulomb\n");
        break;
      }
    }
  else {
    use_ice=0;
    fprintf(echo," No ice\n\n");
    }

  drifters = (tugo_cfg->drifters2D->OnOffFlag.actual_value == KEY_TRUE) ;
  fprintf(echo, "2D drifters simulation %d (flag)\n",drifters);

  tracers = (tugo_cfg->tracers2D->OnOffFlag.actual_value == KEY_TRUE) ;
  fprintf(echo, "2D tracers simulation %d (flag)\n",tracers);

  FlushTracking = (tugo_cfg->misc->FlushTracking.actual_value == KEY_TRUE) ;
  if(FlushTracking) {
    fprintf(echo, " Flush tracking file\n");
    }
  else
    fprintf(echo, " Do NOT flush tracking file\n");

  shear_drag = (tugo_cfg->dissipation->ShearDragFlag.actual_value == KEY_TRUE) ;
  if(shear_drag) {
    fprintf(echo, " Shear drag included in computation\n");
    }
  else
    fprintf(echo, " NO shear drag included in computation\n");

  fprintf(echo,"\n");
/**----------------------------------------------------------------------
  initialization routine, flushing ok*/
  fflush(echo);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void complete_configuration02(tugo_cfg_C *tugo_cfg, int create)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace Keywords ;

  int i, j, n = 499, Days, Hours, Minutes, idum, ndum, hour, nitems;
  int status;

  size_t l;
  char msg[1024];
  char key[500];
  char dum[500];
  char *s, *turbulence_type_string;
  const char *CommentsEnd = "EndComments";
  float Seconds, fdum1, dummy, fdum;
  int ErrorNumber;
  string face_value;

  if(create) create_OuputDirectory(tugo_cfg);

  sprintf(trackingFile, "%s/tracking.dat",gOutputPath);

  fprintf(echo, "Computation mode  : %s \n", tugo_cfg->model->mode.actual_value.c_str()) ;

  fprintf(echo, "Mass conservation : %s \n", tugo_cfg->model->compressibility.face_value().c_str()) ;

  //status = get_char_value(tugo_cfg->dissipation->BottomFrictionType.actual_value, tugo_cfg->dissipation->BottomFrictionType.actual_value) ;
  fprintf(echo, "Bottom friction type: %s \n", tugo_cfg->dissipation->BottomFrictionType.actual_value.c_str()) ;

  if(tugo_cfg->dynamics3D->TimeStep.face_value()=="$2D_time_step") {
    baroclinic_dT=deltaT[0];
    }
  else {
    baroclinic_dT=tugo_cfg->dynamics3D->TimeStep.numerical_value<double>();
    }
  fprintf(echo, " Baroclinic time step:- %f\n", baroclinic_dT);

  if (tugo_cfg->model->compressibility.face_value() == "COMPRESSIBLE") {
    compressibility = COMPRESSIBLE;
    }
  else if (tugo_cfg->model->compressibility.face_value() == "BOUSSINESQ") {
    compressibility = BOUSSINESQ;
    }
  else {
    check_error(-1, "Invalid compressibility mode", __LINE__, __FILE__, 1) ;
    }

/**----------------------------------------------------------------------------
  parse sampling of automatic restart file production*/
  if(tugo_cfg->model->AutomaticInterval.face_value()=="YEARLY") {
    AutomaticInterval=YEARLY;
    }
  else if(tugo_cfg->model->AutomaticInterval.face_value()=="MONTHLY") {
    AutomaticInterval=MONTHLY;
    }
  else if(tugo_cfg->model->AutomaticInterval.face_value()=="WEEKLY") {
    AutomaticInterval=WEEKLY;
    }
  else if(tugo_cfg->model->AutomaticInterval.face_value()=="DAILY") {
    AutomaticInterval=DAILY;
    }
  else if(tugo_cfg->model->AutomaticInterval.face_value()=="HOURLY") {
    AutomaticInterval=HOURLY;
    }
  else {
    check_error(-1, "Invalid automatic interval option", __LINE__, __FILE__, 1) ;
    }

/**----------------------------------------------------------------------------
  parse bottom friction parameterisation*/
  if (tugo_cfg->dissipation->BottomFrictionType.face_value() == KEY_LINEAR) {
    friction_mode=LINEAR_FRICTION;
    }
  else if (tugo_cfg->dissipation->BottomFrictionType.face_value() == KEY_QUADRATIC) {
    friction_mode=QUADRATIC_FRICTION;
    }
  else if (tugo_cfg->dissipation->BottomFrictionType.face_value() == KEY_KARMAN) {
    friction_mode=KARMAN_FRICTION;
    }
  else {
    check_error(-1, "Invalid bottom friction type", __LINE__, __FILE__, 1) ;
    }

  minKv = tugo_cfg->dynamics3D->VerticalDiffusion.numerical_value<double>() ;
  fprintf(echo, "Minimum or constant vertical momentum diffusion coefficient:- %lf\n", minKv);

  minKl = tugo_cfg->dynamics3D->LateralDiffusion.numerical_value<double>() ;
  fprintf(echo, "Minimum or constant lateral momentum diffusion coefficient:- %lf\n", minKl);

  tracers3D_minKh = tugo_cfg->tracers3D->HorizontalDiffusion.numerical_value<double>() ;
  fprintf(echo, "Minimum or constant horizontal tracers diffusion coefficient:- %lf\n", tracers3D_minKh);

  tracers3D_minKv = tugo_cfg->tracers3D->VerticalDiffusion.numerical_value<double>() ;
  fprintf(echo, "Minimum or constant vertical tracers diffusion coefficient:- %lf\n", tracers3D_minKv);

/**----------------------------------------------------------------------------
  parse 3D momentum diffusion oefficient mode*/
  strcpy(dum, tugo_cfg->dissipation->HorizontalViscMode.face_value().c_str()) ;
  if(strcmp(dum,"CONSTANT")==0) {
    Kh_mode=KH_CONSTANT;
    }
  else if(strcmp(dum,"PROPORTIONAL")==0) {
    Kh_mode=KH_PROPORTIONAL;
    }
  else if(strcmp(dum,"SMAGORINSKY")==0) {
    Kh_mode=KH_SMAGORINSKY;
    }
  else if(strcmp(dum,"SMAGORINSKY-2D")==0) {
    Kh_mode=KH_SMAGORINSKY_2D;
    }
  else if(strcmp(dum,"NONE")==0) {
    Kh_mode=KH_NONE;
    }
  else {
    Kh_mode=-1;
    }
  fprintf(echo, "Kh 3D mode: % d (%s) \n",Kh_mode,dum);

/**----------------------------------------------------------------------------
  parse 2D momentum diffusion oefficient mode*/
  strcpy(dum, tugo_cfg->dynamics2D->HorizontalViscMode.face_value().c_str()) ;
  if(strcmp(dum,"CONSTANT")==0) {
    Kh2D_mode=KH_CONSTANT;
    }
  else if(strcmp(dum,"PROPORTIONAL")==0) {
    Kh2D_mode=KH_PROPORTIONAL;
    }
  else if(strcmp(dum,"SMAGORINSKY")==0) {
    Kh2D_mode=KH_SMAGORINSKY;
    }
  else if(strcmp(dum,"UPWINDING")==0) {
    Kh2D_mode=KH_UPWINDING;
    }
  else if(strcmp(dum,"NONE")==0) {
    Kh2D_mode=KH_NONE;
    }
  else {
    Kh2D_mode=-1;
    }
  fprintf(echo, "Kh 2D mode: % d (%s) \n",Kh2D_mode,dum);

  momentum2D_UpwindingCoefficient = tugo_cfg->dynamics2D->UpwindCoefficient.numerical_value<double>() ;
  fprintf(echo, "2D momentum upwinding diffusion coefficient:- %lf\n", momentum2D_UpwindingCoefficient);

/**----------------------------------------------------------------------------
  parse 2D tracer diffusion oefficient mode*/
  strcpy(dum, tugo_cfg->tracers2D->HorizontalDiffusionMode.face_value().c_str()) ;
  if(strcmp(dum,"DEFAULT")==0) {
    Kh2Dtracer_mode=Kh2D_mode;
    }
  else if(strcmp(dum,"NONE")==0) {
    Kh2Dtracer_mode=KH_CONSTANT;
    }
  else if(strcmp(dum,"CONSTANT")==0) {
    Kh2Dtracer_mode=KH_CONSTANT;
    }
  else if(strcmp(dum,"PROPORTIONAL")==0) {
    Kh2Dtracer_mode=KH_PROPORTIONAL;
    }
  else if(strcmp(dum,"SMAGORINSKY")==0) {
    Kh2Dtracer_mode=KH_SMAGORINSKY;
    }
  else if(strcmp(dum,"UPWINDING")==0) {
    Kh2Dtracer_mode=KH_UPWINDING;
    }
  else {
    Kh2Dtracer_mode=-1;
    }
  fprintf(echo, "Kh 2D (tracer) mode: % d (%s) \n",Kh2Dtracer_mode,dum);

  tracers2D_minKh = tugo_cfg->tracers2D->HorizontalDiffusionCoefficient.numerical_value<double>() ;
  fprintf(echo, "Minimum or constant horizontal 2D tracers diffusion coefficient:- %lf\n", tracers2D_minKh);

  tracers2D_UpwindingCoefficient = tugo_cfg->tracers2D->UpwindCoefficient.numerical_value<double>() ;
  fprintf(echo, "2D tracers upwinding diffusion coefficient:- %lf\n", tracers2D_UpwindingCoefficient);

/**----------------------------------------------------------------------------
  parse 3D turbulence closure mode*/

  turbulence_type_string=strdup(tugo_cfg->dynamics3D->TurbulenceClosure.face_value().c_str());
  printf("Turbulence type id: %d %s \n",__LINE__,__FILE__);
  if(strcmp(turbulence_type_string,"CONSTANT")==0){
          //tugo_cfg->dynamics3D->TurbulenceClosure.face_value() == KV_STRG_CONSTANT) {
    turbulence_type=KV_CONSTANT;
    }
  else if(strcmp(turbulence_type_string,"GASPAR")==0){
          //if(tugo_cfg->dynamics3D->TurbulenceClosure.face_value() == KV_STRG_GASPAR) {
    turbulence_type=KV_GASPAR;
    }
  else if(strcmp(turbulence_type_string,"HOMOGENEOUS")==0){
          //if(tugo_cfg->dynamics3D->TurbulenceClosure.face_value() == KV_STRG_HOMOGENEOUS) {
    turbulence_type=KV_HOMOGENEOUS;
    }
  else if(strcmp(turbulence_type_string,"TKE")==0){
                  //if(tugo_cfg->dynamics3D->TurbulenceClosure.face_value() == KV_STRG_TKE) {
    turbulence_type=KV_TKE;
    printf("Turbulence type id: %d  \n",turbulence_type);
    }
  else {
    check_error(-1, "Invalid turbulence closure mode", __LINE__, __FILE__, 1) ;
    }
  fprintf(echo, "Turbulence type: %d (%s) \n",turbulence_type,turbulence_type_string);

  gAsselinC1 = tugo_cfg->dynamics2D->AsselinCoefficientC1.numerical_value<double>() ;
  gAsselinC2 = 1-2*gAsselinC1 ;
  fprintf(echo, "Asselin filter coefficient: %lf (on u0 and u2) %lf (on u1) \n", gAsselinC1, gAsselinC2);

  equilibrium = (tugo_cfg->tides->equilibrium.actual_value == KEY_TRUE);
  if(equilibrium) {
    fprintf(echo, " Equilibrium extension allowed for tidal spectrum \n");
  } else
    fprintf(echo, " NO equilibrium extension allowed for tidal spectrum\n");


/**----------------------------------------------------------------------------
  imbricated boundary conditions*/
  strcpy(dum, tugo_cfg->boundaries->ArchivedBCTransport.face_value().c_str()) ;
  if(strcmp(dum,"u")==0) {
    external_transport=0;
    fprintf(echo, "use currents from the external simulation\n");
    }
  else if(strcmp(dum,"U")==0) {
    external_transport=1;
    fprintf(echo, "use transport from the external simulation\n");
    }
  else if(strcmp(dum,"NULL")==0) {
    external_transport=2;
    fprintf(echo, "use no currents and no transport from the external simulation\n");
    }
  else {
    external_transport=-1;
    }

/**----------------------------------------------------------------------------
  Ensure that background velocity is not null */
  if (gU0 <= 0) gU0 = u0_min ;

/**----------------------------------------------------------------------------
  Analyse the advection-diffusion mode, obsolete, should be handle by g2D_discretisation*/
  if(tugo_cfg->dynamics2D->AdvectionDiffusionDiscretisation.face_value()      == "NODE-WISE") {
    AD_discretisation=NODE_WISE;
    }
  else if(tugo_cfg->dynamics2D->AdvectionDiffusionDiscretisation.face_value() == "ELEMENT-WISE") {
    AD_discretisation=ELEMENT_WISE;
    }
  else {
    AD_discretisation=-1;
    }

/**----------------------------------------------------------------------------
  Analyse the 2D velocity discretisation*/
  face_value=tugo_cfg->dynamics2D->V2DDiscretisation.face_value();
  if(face_value      == KEY_LGP0) {
    u2D_discretisation=LGP0;
    }
  else if(face_value == KEY_LGP1) {
    u2D_discretisation=LGP1;
    }
  else if(face_value == KEY_LGP1_CG) { //for hand-written files
    u2D_discretisation=LGP1;
    }
  else if(face_value == KEY_NCP1) {
    u2D_discretisation=NCP1;
    }
  else if(face_value == KEY_NCP1_CG) {
    u2D_discretisation=NCP1;
    }
  else if(face_value == KEY_LGP1_DG) {
    u2D_discretisation=DGP1;
    }
  else if(face_value == KEY_NCP1_DG) {
    u2D_discretisation=DNP1;
    }
  else if(face_value == KEY_CQP0) {
    u2D_discretisation=CQP0;
    }
  else if(face_value == KEY_CQN1) {
    u2D_discretisation=CQN1;
    }
  else if(face_value == KEY_CQP1) {
    u2D_discretisation=CQP1;
    }
  else {
    u2D_discretisation=-1;
    }

/**----------------------------------------------------------------------------
  Analyse the 2D pressure (elevation) discretisation*/
  face_value=tugo_cfg->dynamics2D->SLVDiscretisation.face_value();
  if(face_value      == KEY_LGP1) {
    z2D_discretisation=LGP1;
    }
  else if(face_value == KEY_LGP1_CG) {
    z2D_discretisation=LGP1;
    }
  else if(face_value == KEY_LGP2_CG) {
    z2D_discretisation=LGP2;
    }
  else if(face_value == KEY_LGP0) {
    z2D_discretisation=LGP0;
    }
  else if(face_value == KEY_CQP0) {
    z2D_discretisation=CQP0;
    }
  else if(face_value == KEY_CQP1) {
    z2D_discretisation=CQP1;
    }
  else {
    z2D_discretisation=-1;
    }

/**----------------------------------------------------------------------------
  Analyse the 2D velocity gradient discretisation*/
  if(tugo_cfg->dynamics2D->G2DDiscretisation.face_value()      == KEY_LGP0) {
    g2D_discretisation=LGP0;
    }
  else if(tugo_cfg->dynamics2D->G2DDiscretisation.face_value() == KEY_LGP1) {
    g2D_discretisation=LGP1;
    }
  else if(tugo_cfg->dynamics2D->G2DDiscretisation.face_value() == KEY_NCP1) {
    g2D_discretisation=NCP1;
    }
  else if(tugo_cfg->dynamics2D->G2DDiscretisation.face_value() == KEY_NCP1_CG) {
    g2D_discretisation=NCP1;
    }
  else if(tugo_cfg->dynamics2D->G2DDiscretisation.face_value() == KEY_LGP1_DG) {
    g2D_discretisation=DGP1;
    }
  else if(tugo_cfg->dynamics2D->G2DDiscretisation.face_value() == KEY_NCP1_DG) {
    g2D_discretisation=DNP1;
    }
  else if(tugo_cfg->dynamics2D->G2DDiscretisation.face_value() == KEY_LGP1) {
    g2D_discretisation=LGP2;
    }
  else if(tugo_cfg->dynamics2D->G2DDiscretisation.face_value() == KEY_INTRINSIC) {
    g2D_discretisation=INTRINSIC;
    }
  else if(tugo_cfg->dynamics2D->G2DDiscretisation.face_value() == KEY_AUTO) {
    g2D_discretisation=AUTO;
    }
  else {
    g2D_discretisation=-1;
    }

/**----------------------------------------------------------------------------
  at the right place?*/
  if(g2D_discretisation==AUTO) {
    g2D_discretisation=gradient_natural_discretisation(u2D_discretisation);
    }
    
/**----------------------------------------------------------------------------
  transition patch*/
  switch(AD_discretisation) {
    case ELEMENT_WISE:
//       if(g2D_discretisation!=INTRINSIC) {
//         char *msg;
//         asprintf(&msg,"WARNING : u-gradient discretisation is not set to INTRINSIC, may conflict with advection/diffusion mode setting");
//         check_error(1, msg, __LINE__, __FILE__, 0);
//         free(msg);
//         }
//       g2D_discretisation=LGP0;
      break;
    case NODE_WISE:
      if(g2D_discretisation!=INTRINSIC) {
        char *msg;
        asprintf(&msg,"WARNING : u-gradient discretisation is not set to INTRINSIC, may conflict with advection/diffusion mode setting");
        check_error(1, msg, __LINE__, __FILE__, 0);
        free(msg);
        }
      g2D_discretisation=u2D_discretisation;
      break;
    }
    
/**----------------------------------------------------------------------------
  mysteriously disappeared after 2011-05-04, re-intalled 2012-07-06 */
  if (tugo_cfg->dynamics2D->Advection2Dmode.face_value()      == "CONSERVATIVE") {
    advection2D_mode=CONSERVATIVE;
    }
  else if (tugo_cfg->dynamics2D->Advection2Dmode.face_value() == "NON-CONSERVATIVE") {
    advection2D_mode=PRIMITIVE;
    }
  else if (tugo_cfg->dynamics2D->Advection2Dmode.face_value() == "LAGRANGIAN") {
    advection2D_mode=LAGRANGIAN;
    }
  else {
    advection2D_mode=-1;
    }

/**----------------------------------------------------------------------------
  Analyse the 3D velocity discretisation*/
  switch (u2D_discretisation) {
    case LGP0:
      tugo_cfg->dynamics3D->VDiscretisation.actual_value = (string) "LGP0";
      break;
    case LGP1:
      tugo_cfg->dynamics3D->VDiscretisation.actual_value = (string) "LGP1";
      break;
    case NCP1:
      tugo_cfg->dynamics3D->VDiscretisation.actual_value = (string) "NCP1";
      break;
    }

  gSequentialPaire2D=paire_definition(z2D_discretisation,u2D_discretisation);

  if(tugo_cfg->dynamics2D->WaveScalarProduct.face_value() == KEY_NQUAD) {
    sproduct2D=NQUAD;
    }
  else {
    sproduct2D=INTGL;
    }

  if (tugo_cfg->dynamics2D->Momentum2DMode.face_value()      == "VELOCITY") {
    momentum2d_mode=VELOCITY;
    }
  else if (tugo_cfg->dynamics2D->Momentum2DMode.face_value() == "TRANSPORT") {
    momentum2d_mode=TRANSPORT;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void check_configuration(tugo_cfg_C *tugo_cfg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace Keywords ;

  int i, j, n = 499, Days, Hours, Minutes, idum, ndum, hour, nitems;
  int status;

  size_t l;
  char msg[1024];
  char key[500];
  char dum[500];
  char *s;
  const char *CommentsEnd = "EndComments";
  float Seconds, fdum1, dummy, fdum;

  switch (we_type) {
    case (GWE_CLX):
/*----------------------------------------------------------------------
      GWE, nodal quadrature matrix LGP1*/
      tugo_cfg->dynamics2D->Continuity2DIntegration.actual_value = "LEAPFROG";
      tugo_cfg->dynamics2D->Momentum2DIntegration.actual_value   = "FORWARD";
      break;

    case (GWE_XXX):
/*----------------------------------------------------------------------
      GWE new formulation, nodal quadrature matrix LGP1*/
      break;

    case (WE_EXPLICIT):
/*----------------------------------------------------------------------
      WE explicit,*/
      if(gTheta != 0.) {
//        check_error(-1, "theta should be zero, but is not in input file", __LINE__, __FILE__, 1);
        }
      break;

    case (WE_IMPLICIT):
/*----------------------------------------------------------------------
      WE implicit*/
      break;

    case (FV_IMPLICIT):
/*----------------------------------------------------------------------
      Finite Volumes implicit*/
      break;

    default:
      check_error(-1, "illicit shallow-water mode", __LINE__, __FILE__, 1);
      break;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int configure_init(char *filename, int fmt, tugo_cfg_C *tugo_cfg, int delta)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  string out;
  int create=!delta;

  if(fmt==0) {
    /* STS: NO old input format */
    }
  else {
    hugo_readinput(tugo_cfg,filename);
    /* STS: keep gOutputPath */
    if(!gOutputPath)
      gOutputPath=strdup(tugo_cfg->model->OutputPath.face_value().c_str());
    complete_configuration01(tugo_cfg, create);
    }

  complete_configuration02(tugo_cfg,create);
  check_configuration(tugo_cfg);

  return(0);
}
