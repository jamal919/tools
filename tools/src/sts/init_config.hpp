
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

//  init_config.hpp 
//  
//  
//  Toutes ces CLASSES sont utiles pour lire et ecrire les fichiers de configuration
//  du modele TUGO
//  
//  La classe la plus basses est KEY_C
//  et la classe la plus haute tugo_cfg_C 
//
//  le fichier est indissociable de init_config_constructor.cpp
//
//  Thierry LETELLIER LEGOS Toulouse 08Jun2006
//

#ifndef _UGO_CONFIGURE_H
#define _UGO_CONFIGURE_H
#include <iostream>
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include <cmath>
#include <string.h> 
#include <string> 
#include "exceptions.hpp"
#include "keywords.hpp"
#include <functions.h>
#include <poc-array.hpp>

using namespace std;

class KEY_C {
public :
  char *keyword;
  char *default_value;
  //char *actual_value;
  string actual_value ;

  char *realized_value;
  bool imported;
  bool initialised;
  char *comments;
  char *more_comments;
  char *choice;
  int NBchoice;
  int type;
  char *units;
  bool deprecated;

  KEY_C() {
     keyword        = new char[128];
     default_value  = new char[128];

//     actual_value   = new char[128];

     realized_value = 0;
     imported       = false;
     initialised    = false;
     comments       = new char[512];
     more_comments  = 0;
     choice         = new char[128];
     NBchoice       = 0;
     units          = new char[56] ;
     type           = -1;
     deprecated     = false;
     }

   int destruct () {
     deletep(&keyword);
     deletep(&default_value);
     deletep(&realized_value);
     deletep(&comments);
     deletep(&more_comments);
     deletep(&choice);
     deletep(&units);
     }

  void  implement(const char *A, const char *B, const char *C, const char *D, int E , const char *F, const char *G){
    sprintf(keyword,"%s",A);
    sprintf(default_value,"%s",B);
    if(B!=C) {
      printf("warning, keyword %s : default=%s but initialisation=%s\n",A,B,C);
      }
    actual_value = C;

    sprintf(comments,"%s",D);
    NBchoice=E;
    sprintf(choice,"%s",F);
    sprintf(units,"%s",G);
    }

  void  implement_more_comments(const char *A){
    int longeur;
    longeur = strlen(A)+1;
    more_comments = new char[longeur];
    sprintf(more_comments,"%s",A);
    }

  // --- Transmit actual_value but without "<default>".

  string face_value (void) const throw () { 
    string s (actual_value) ;
    if(deprecated) {
      cout <<__FUNCTION__<< ", attempt to use deprecated parameter " << keyword << "\n" << flush;
      exit(-1);
      }
    // skip <default> if present
    const int pos = s.find("<default>") ;
    if (pos >= 0) s = s.substr(pos+9) ;

    // j'enleve les espaces a la fin
    //while (s.rfind(" ") == s.size() - 1) s = s.erase(s.size() - 1, 1) ;
    /* The previous line seems not to work : everything append as if s.size()
       in the condition is not re-evaluate during the loop. */
    bool cont = true ;
    while (cont) {
      const string::size_type ending = s.size() - 1 ;
      if (s.rfind(" ") == ending) s.erase(ending, 1) ;
      else cont = false ;
      }

    // enfin j'enleve les espaces au debut
    while (s.find(" ") == 0) s = s.erase(0, 1) ;

    return s ;
    }   // string face_value (void)

  /* --- Overloaded function, which transmits the numerical value of
         actual_value, if it is a number stocked as a string, or an exception. */

  template <class T> inline T numerical_value (void) const throw () {
    if(face_value().c_str()==Keywords::KEY_TRUE)  return 1;
    if(face_value().c_str()==Keywords::KEY_FALSE) return 0;
    string s=face_value();
    T x=atof(s.c_str());
    return static_cast<T>(x) ;
  }   // T numerical_value (void)

  template <class T> inline T multiple_value (int pos) const throw () {
    string s (face_value().c_str());
    char *actual,*bkp;
    char *next;
    int k=0;
    T value;
    actual=strdup(face_value().c_str());
    bkp=actual;
    while(k<pos) {
      strtod(actual,&next);
      k++;
      actual=next;
      }
    value=strtod(actual,&next);
    free(bkp);
    return static_cast<T>(value) ;
  }   // T multiple_value (int)

  inline int numerical_value (void) const throw () {
    string s=face_value();
    int x=atoi(s.c_str());
    return x ;
  }   // int numerical_value (void)

private :
  // no private members 

};


//##################################
//
// class of model Key Word ....
//
//##################################
class model_C {
public :
  KEY_C rootname;
  KEY_C TimeStep, SubTimeStep;
  KEY_C OrigineTime;
  KEY_C MeshType,MeshFormat,MeshMeta;
  KEY_C MeshFile;
  
  KEY_C TopoFile,SlopeFile,ZminFile;
  KEY_C MinDepth,RelativeMinDepth;
  
  KEY_C BelFile;
  KEY_C RunDuration;
  KEY_C SpinUpDuration, SpinUpMode;
  KEY_C RestartFile, RestartFormat;
  KEY_C ContinuationFile, ContinuationFormat, ContinuationInterval, AutomaticInterval;
  KEY_C OutputPath;
  KEY_C LinearSolver, SubLinearSolver, SolverMode;
  KEY_C mode;
  KEY_C compressibility;
  KEY_C WaveEqForm;
  KEY_C TwoDMometumEqForm;
  KEY_C BarotropicCurrentsDiscretisation;
  KEY_C tau0, theta;
  KEY_C OPENMP_nCPUsMax;
  KEY_C MetisPath;
  KEY_C ThreeDPressure;
  KEY_C EchoFile;
  KEY_C CompatibilityScale;

  //constructor puts default value of  variables
  model_C();
  int destruct () {
   rootname.destruct();
   TimeStep.destruct();
   SubTimeStep.destruct();
   OrigineTime.destruct();
   MeshType.destruct();
   MeshFormat.destruct();
   MeshMeta.destruct();
   MeshFile.destruct();
   
   TopoFile.destruct();
   SlopeFile.destruct();
   ZminFile.destruct();
   MinDepth.destruct();
   RelativeMinDepth.destruct();
   
   BelFile.destruct();
   RunDuration.destruct();
   SpinUpDuration.destruct();
   SpinUpMode.destruct();
   RestartFile.destruct();
   RestartFormat.destruct();
   ContinuationFile.destruct();
   ContinuationFormat.destruct();
   ContinuationInterval.destruct();
   AutomaticInterval.destruct();
   OutputPath.destruct();
   LinearSolver.destruct();
   SubLinearSolver.destruct();
   SolverMode.destruct();
   mode.destruct();
   compressibility.destruct();
   WaveEqForm.destruct();
   TwoDMometumEqForm.destruct();
   BarotropicCurrentsDiscretisation.destruct();
   tau0.destruct();
   theta.destruct();
   OPENMP_nCPUsMax.destruct();
   MetisPath.destruct();
   ThreeDPressure.destruct();
   EchoFile.destruct();
   CompatibilityScale.destruct();
    }

private :
 // no private members 
};

//##################################
//
// class of physics Key Word ....
//
//##################################
class physics_C {
public :
  KEY_C gravity_mode, anelastic_body_tide;
  KEY_C gravity_adjustment;
  KEY_C coriolis_mode, coriolis_latitude;
  KEY_C water_viscosity;

  //constructor puts default value of  variables
  physics_C();
  int destruct () {
    gravity_mode.destruct();
    anelastic_body_tide.destruct();
    gravity_adjustment.destruct();
    coriolis_mode.destruct();
    coriolis_latitude.destruct();
    water_viscosity.destruct();
    }

private :
 // no private members 
};

//##################################
//
// class of dissipation Key Word ....
//
//##################################
class dissipation_C {
public :
  KEY_C MinBackGroundSpd;
  KEY_C BackgroundVelocityFile;
  KEY_C BKGPolygons,BKGValues;
  KEY_C BottomFrictionType;
  KEY_C BottomRugosityFile;
  KEY_C IceShelvesFile,IceShelvesCoefficient;
  KEY_C LinearFrictionCoeff;
  KEY_C QuadraticFrictionCoeff;
  KEY_C BottomFrictionCoeffFile;
  KEY_C QFCPolygons,QFCValues;
  KEY_C KarmanRugosity,RoughnessLength;
  KEY_C MinQuadraticFrictionCoeff;
  KEY_C QFCMPolygons,QFCMValues;
  KEY_C RugosityPolygons,RugosityValues;
  KEY_C RugosityValuesbyRegions;
  KEY_C InternDragAlgorithm;
  KEY_C InternDragSlope;
  KEY_C InternDragHmin;
  KEY_C InternDragHmax;
  KEY_C IWD01Polygons,IWD01Values;
  KEY_C IWD01ValuesByRegions;
  KEY_C InternDragRugosity;
  KEY_C MixedLayerCoeff;
  KEY_C SmagorinskyCoefficient,MinhorizontalVisc,HorizontalViscMode;
  KEY_C ShearDragFlag;
  KEY_C BruntVassalaValue;
  KEY_C BruntVassalaFile;

  //constructor puts default value of  variables
  dissipation_C();

  int destruct () {
   MinBackGroundSpd.destruct();
   BackgroundVelocityFile.destruct();
   BKGPolygons.destruct();
   BKGValues.destruct();
   BottomFrictionType.destruct();
   IceShelvesFile.destruct();
   IceShelvesCoefficient.destruct();
   BottomRugosityFile.destruct();
   LinearFrictionCoeff.destruct();
   QuadraticFrictionCoeff.destruct();
   BottomFrictionCoeffFile.destruct();
   QFCPolygons.destruct();
   QFCValues.destruct();
   KarmanRugosity.destruct();
   RoughnessLength.destruct();
   MinQuadraticFrictionCoeff.destruct();
   QFCMPolygons.destruct();
   QFCMValues.destruct();
   RugosityPolygons.destruct();
   RugosityValues.destruct();
   RugosityValuesbyRegions.destruct();
   InternDragAlgorithm.destruct();
   InternDragSlope.destruct();
   InternDragHmin.destruct();
   InternDragHmax.destruct();
   IWD01Polygons.destruct();
   IWD01Values.destruct();
   IWD01ValuesByRegions.destruct();
   InternDragRugosity.destruct();
   MixedLayerCoeff.destruct();
   SmagorinskyCoefficient.destruct();
   MinhorizontalVisc.destruct();
   HorizontalViscMode.destruct();
   ShearDragFlag.destruct();
   BruntVassalaValue.destruct();
   BruntVassalaFile.destruct();
  }
private :
 // no private members 
};

//##################################
//
// class of tide Key Word ....
//
//##################################
class tides_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C AstronomicPotentialFlag;
  KEY_C LSAFlag;
  KEY_C LSADir;
  KEY_C LSAConvention;
  KEY_C AtlasDir;
  KEY_C AtlasConvention;
  KEY_C PressureFlag;
  KEY_C PressureDir;
  KEY_C PressureConvention;
  KEY_C admittance;
  KEY_C equilibrium;
  KEY_C NodalCorrection;
  KEY_C LSACoeff;
  KEY_C DeformationLoveNum;
  KEY_C BoundTideFlag;
  KEY_C BoundTideFile;
  KEY_C Spectral2DRun, SpectralSolver;
  KEY_C DominantWave_1, DominantWave_2, PriorSolution;
  KEY_C SpectralPaire;
  KEY_C Spectral2DMaxIteration, CompoundMaxIteration;
  KEY_C SpecificData, TopoFile, SlopeFile;
  KEY_C SpectralRecycle;
  KEY_C ReducedArchive;
  KEY_C SpectralArchiveLevel;
  KEY_C Spectral3DRun, Spectral3DOBCs;
  KEY_C VDiffusionTimeStep, VDiffusionDuration, Spectral3DMaxIteration, Spectral3DFirstIteration;
  bool  active;

  tides_C();

  int destruct () {
   OnOffFlag.destruct();
   AstronomicPotentialFlag.destruct();
   LSAFlag.destruct();
   LSADir.destruct();
   LSAConvention.destruct();
   AtlasDir.destruct();
   AtlasConvention.destruct();
   PressureFlag.destruct();
   PressureDir.destruct();
   PressureConvention.destruct();
   admittance.destruct();
   equilibrium.destruct();
   NodalCorrection.destruct();
   LSACoeff.destruct();
   DeformationLoveNum.destruct();
   BoundTideFlag.destruct();
   BoundTideFile.destruct();
   Spectral2DRun.destruct();
   SpectralSolver.destruct();
   DominantWave_1.destruct();
   DominantWave_2.destruct();
   PriorSolution.destruct();
   SpectralPaire.destruct();
   Spectral2DMaxIteration.destruct();
   CompoundMaxIteration.destruct();
   SpecificData.destruct();
   TopoFile.destruct();
   SlopeFile.destruct();
   SlopeFile.destruct();
   SpectralRecycle.destruct();
   ReducedArchive.destruct();
   SpectralArchiveLevel.destruct();
   Spectral3DRun.destruct();
   Spectral3DOBCs.destruct();
   VDiffusionDuration.destruct();
   VDiffusionTimeStep.destruct();
   Spectral3DFirstIteration.destruct();
   Spectral3DMaxIteration.destruct();
  }

private :
 // no private members 
};
//##################################
//
// class of atmosphere Key Word ....
//
//##################################
class atmosphere_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C PressureFlag;
  KEY_C PressureDir;
  KEY_C PressureConvention;
  KEY_C WindFlag;
  KEY_C LongWave,ShortWave;
  KEY_C AtmosDirectory;
  KEY_C FileFormat;
  KEY_C FileEndian;
  KEY_C FileConvention;
  KEY_C ApplyMask, MaskFile;
  KEY_C KeepMeanWindFlag;
  KEY_C KeepMeanPressureFlag;
  KEY_C TimeInterp;
  KEY_C BulkFormula,WindStress;
  KEY_C BoundInversebarometer;
  bool  active;

  atmosphere_C();

  int destruct () {
   OnOffFlag.destruct();
   PressureFlag.destruct();
   PressureFlag.destruct();
   PressureDir.destruct();
   PressureConvention.destruct();
   WindFlag.destruct();
   LongWave.destruct();
   ShortWave.destruct();
   AtmosDirectory.destruct();
   FileFormat.destruct();
   FileEndian.destruct();
   FileConvention.destruct();
   ApplyMask.destruct();
   MaskFile.destruct();
   KeepMeanWindFlag.destruct();
   KeepMeanPressureFlag.destruct();
   TimeInterp.destruct();
   BulkFormula.destruct();
   WindStress.destruct();
   BoundInversebarometer.destruct();
  }

private :
 // no private members 
};

//##################################
//
// class of ocean wave coupling Key Word ....
//
//##################################
class OceanWaves_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C Parameterisation;
  KEY_C DragFlag;
  KEY_C FrictionFlag;
  KEY_C TimeInterp;
  KEY_C Discretisation;
  KEY_C FileDirectory;
  KEY_C FileFormat;
  KEY_C FileConvention;
  KEY_C MaskFile;
  bool  active;

  OceanWaves_C();
  
  int destruct () {
   OnOffFlag.destruct();
   Parameterisation.destruct();
   DragFlag.destruct();
   FrictionFlag.destruct();
   TimeInterp.destruct();
   Discretisation.destruct();
   FileDirectory.destruct();
   FileFormat.destruct();
   FileConvention.destruct();
   MaskFile.destruct();
  }
private :
 // no private members 
};

//##################################
//
// class of boundaries Key Word ....
//
//##################################
class boundaries_C {
public :
  //constructor puts default value of  variables
  KEY_C type,mode;
  KEY_C RelaxTime;
  KEY_C ArchivedBCUse;
  KEY_C ArchivedBCDirectory;
  KEY_C ArchivedBCFormat;
  KEY_C ArchivedBCRoot;
  KEY_C ArchivedBCUnits;
  KEY_C ArchivedBCTransport;
  KEY_C PeriodicFile;
  KEY_C SpongeFlag;
  KEY_C PMLfactor,BufferZoneFile;
  KEY_C RiversFlag,RiversFile,RiversTimeTemplate;
  KEY_C SolidCondition;
  KEY_C NoSlipSize;
  bool  active;

  boundaries_C();

  int destruct () {
    type.destruct();
    mode.destruct();
    RelaxTime.destruct();
    
    ArchivedBCUse.destruct();
    ArchivedBCDirectory.destruct();
    ArchivedBCFormat.destruct();
    ArchivedBCRoot.destruct();
    ArchivedBCUnits.destruct();
    ArchivedBCTransport.destruct();
    
    PeriodicFile.destruct();
    SpongeFlag.destruct();
    PMLfactor.destruct();
    BufferZoneFile.destruct();
    
    RiversFlag.destruct();
    RiversFile.destruct();
    RiversTimeTemplate.destruct();
    SolidCondition.destruct();
    NoSlipSize.destruct();
    }

private :
 // no private members 
};

//##################################
//
// class of ice Key Word ....
//
//##################################
class ice_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C IceCoverFile;
  KEY_C IceElasticityFile;
  KEY_C TimeStep;
  KEY_C MomentumSolver;
  KEY_C Rheology;
  bool  active;

  ice_C();

  int destruct () {
    OnOffFlag.destruct();
    IceCoverFile.destruct();
    IceElasticityFile.destruct();
    TimeStep.destruct();
    MomentumSolver.destruct();
    Rheology.destruct();
    }

private :
 // no private members 
};

//##################################
//
// class of topography Key Word ....
//
//##################################
class topography_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C BottomMotionFlag;
  KEY_C BottomMotionDirectory;
  KEY_C BottomMotionFileConvent;
  
  KEY_C TopoFile,SlopeFile,ZminFile;
  KEY_C MinDepth,RelativeMinDepth;
  KEY_C DryingMinDepth;

  KEY_C DepthOffsetPolygons, DepthOffsetValues;
  KEY_C DepthScaleFactor;
  KEY_C MeanLevelFile;
  bool  active;
  bool  BottomMotionActive;

  topography_C();

  int destruct () {
    OnOffFlag.destruct();
    BottomMotionFlag.destruct();
    BottomMotionDirectory.destruct();
    BottomMotionFileConvent.destruct();
    
    TopoFile.destruct();
    SlopeFile.destruct();
    ZminFile.destruct();
    MinDepth.destruct();
    RelativeMinDepth.destruct();
    DryingMinDepth.destruct();
    
    DepthOffsetPolygons.destruct();
    DepthOffsetValues.destruct();
    DepthScaleFactor.destruct();
    MeanLevelFile.destruct();
    }

private :
 // no private members 
};
//##################################
//
// class of misc Key Word ....
//
//##################################
class misc_C {
public :
  //constructor puts default value of  variables
  KEY_C NonTidalLoadFlag;
  KEY_C SmoothingFlag;
  KEY_C checks;
  KEY_C Check_Wconsistency;
  KEY_C CompatibilityScale;
  KEY_C BarotropicW;
  KEY_C coupling_2D_3D;
  KEY_C StabilityCheck,StabilityControl;
  KEY_C MaxDeltaElevation,MaxDeltaCurrents,MaxDeltaPeriod;
  KEY_C nSubCycles;
  KEY_C SubcyclingPersistence;
  KEY_C SubcyclingTemporization;
  KEY_C SubcyclingExtension;
  KEY_C SubcyclingLockLimit;
  KEY_C SubcyclingCFLratio;
  KEY_C SmoothZField,SmoothUField,SmoothAdvection,SmoothDiffusion;
  KEY_C FlushTracking;
  KEY_C ToponymsFile;
  bool OnOffFlag;

  misc_C();
  int destruct () {
   NonTidalLoadFlag.destruct();
   SmoothingFlag.destruct();
   checks.destruct();
   Check_Wconsistency.destruct();
   CompatibilityScale.destruct();
   BarotropicW.destruct();
   coupling_2D_3D.destruct();
   StabilityCheck.destruct();
   StabilityControl.destruct();
   MaxDeltaElevation.destruct();
   MaxDeltaCurrents.destruct();
   MaxDeltaPeriod.destruct();
   nSubCycles.destruct();
   SubcyclingPersistence.destruct();
   SubcyclingTemporization.destruct();
   SubcyclingExtension.destruct();
   SubcyclingLockLimit.destruct();
   SubcyclingCFLratio.destruct();
   SmoothZField.destruct();
   SmoothUField.destruct();
   SmoothAdvection.destruct();
   SmoothDiffusion.destruct();
   FlushTracking.destruct();
   ToponymsFile.destruct();
   }

private :
 // no private members 
};
//##################################
//
// class of fearch Key Word ...
//
//##################################
class fearchive_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C StateFile, VectorFile;
  KEY_C ArchFormat, ArchConvention;
  KEY_C ArchStart,ArchInterval;
  KEY_C ArchBits;
  KEY_C Arch_huv_Scale,Arch_pW_Scale;
  KEY_C ArchOverAppendFlag;
  KEY_C OverwriteAppendFlag;
  KEY_C z_flag,ubar_flag;
  KEY_C transport2D,divergence2D;
  KEY_C ws_flag,Pa_flag,ibd_flag;
  KEY_C u_flag,w_flag,T_flag,S_flag,R_flag;
  KEY_C horizontal_diffusion2D,horizontal_advection2D;
  KEY_C horizontal_diffusion,horizontal_advection;
  KEY_C Cd_flag;
  KEY_C lateral_mixing;
  KEY_C radiation_stress;
  KEY_C gradSxx;
  bool  active;

  fearchive_C();

  int destruct () {
   OnOffFlag.destruct();
   StateFile.destruct(); VectorFile.destruct();
   ArchFormat.destruct(); ArchConvention.destruct();
   ArchStart.destruct();ArchInterval.destruct();
   ArchBits.destruct();
   Arch_huv_Scale.destruct();Arch_pW_Scale.destruct();
   ArchOverAppendFlag.destruct();
   OverwriteAppendFlag.destruct();
   z_flag.destruct();ubar_flag.destruct();
   transport2D.destruct();divergence2D.destruct();
   ws_flag.destruct();Pa_flag.destruct();ibd_flag.destruct();
   u_flag.destruct();w_flag.destruct();T_flag.destruct();S_flag.destruct();R_flag.destruct();
   horizontal_diffusion2D.destruct();horizontal_advection2D.destruct();
   horizontal_diffusion.destruct();horizontal_advection.destruct();
   Cd_flag.destruct();
   lateral_mixing.destruct();
   radiation_stress.destruct();
   gradSxx.destruct();
    }

private :
 // no private members 
};
//##################################
//
// class of structured Key Word ....
//
//##################################
class structured_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C zone;
  KEY_C ArchiveStart, ArchiveInterval;
  KEY_C Start, Interval;
  KEY_C OverwriteAppendFlag;
  KEY_C z_flag,ubar_flag,Cd_flag,z0_flag,ibd_flag;
  KEY_C MKE2D,zrange;
  KEY_C transport2D,divergence2D;
  KEY_C horizontal_diffusion2D,horizontal_advection2D;
  KEY_C ws_flag,Pa_flag;
  KEY_C u_flag,w_flag,T_flag,S_flag,R_flag;
  KEY_C horizontal_diffusion,pressure_gradient;
  KEY_C horizontal_advection;
  KEY_C lateral_mixing,tracer2D;
  KEY_C radiation_stress;
  bool  active;

  structured_C();
  int destruct () {
   OnOffFlag.destruct();
   zone.destruct();
   ArchiveStart.destruct(); ArchiveInterval.destruct();
   Start.destruct(); Interval.destruct();
   OverwriteAppendFlag.destruct();
   z_flag.destruct();ubar_flag.destruct();Cd_flag.destruct();z0_flag.destruct();ibd_flag.destruct();
   MKE2D.destruct();zrange.destruct();
   transport2D.destruct();divergence2D.destruct();
   horizontal_diffusion2D.destruct();horizontal_advection2D.destruct();
   ws_flag.destruct();Pa_flag.destruct();
   u_flag.destruct();w_flag.destruct();T_flag.destruct();S_flag.destruct();R_flag.destruct();
   horizontal_diffusion.destruct();pressure_gradient.destruct();
   horizontal_advection.destruct();
   lateral_mixing.destruct();tracer2D.destruct();
   radiation_stress.destruct();
     }

private :
 // no private members 
};
//##################################
//
// class of sample Key Word ....
//
//##################################
class sample_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C file;
  KEY_C start;
  KEY_C SaveInterval;
  KEY_C OverAppendFlag;
  bool  active;

  KEY_C SectionOnOffFlag;
  KEY_C SectionInputFile;
  KEY_C SectionStart;
  KEY_C SectionSaveInterval;
  KEY_C SectionOverAppendFlag;
  bool  SectionActive;

  sample_C();
  int destruct () {
   OnOffFlag.destruct();
   file.destruct();
   start.destruct();
   SaveInterval.destruct();
   OverAppendFlag.destruct();
   SectionOnOffFlag.destruct();
   SectionInputFile.destruct();
   SectionStart.destruct();
   SectionSaveInterval.destruct();
   SectionOverAppendFlag.destruct();
    }
private :
 // no private members 
};

//##################################
//
// class of analysis Key Word ....
//
//##################################
class analysis_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C ListFile;
  KEY_C AutoCompletion;
  KEY_C HarmoStart;
  KEY_C SampleInterval;
  KEY_C ComputeInterval;
  KEY_C Strategy;
  KEY_C SequentialRecycle;
  KEY_C SequentialArchiveLevel;
  KEY_C SGArchiveFlag,LGP1ArchiveFlag,NETCDFArchiveFlag;
  KEY_C grid,z_flag,u_flag,LSA_flag;
  bool  active,SGsave;

  analysis_C();
  int destruct () {
    OnOffFlag.destruct();
    ListFile.destruct();
    AutoCompletion.destruct();
    HarmoStart.destruct();
    SampleInterval.destruct();
    ComputeInterval.destruct();
    Strategy.destruct();   
    SequentialRecycle.destruct();
    SequentialArchiveLevel.destruct();
    SGArchiveFlag.destruct();
    LGP1ArchiveFlag.destruct();
    NETCDFArchiveFlag.destruct();
    grid.destruct();z_flag.destruct();u_flag.destruct();LSA_flag.destruct();
    }

private :
 // no private members 
};

//##################################
//
// class of energy Key Word ....
//
//##################################
class energy_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C start;
  KEY_C interval;
  KEY_C TimeStepComp;
  KEY_C HarmoComp, HarmoCompByRegions;
  KEY_C ResonanceExploration;
  KEY_C ResonanceFrequencies;
  bool  active;

  energy_C();
  int destruct () {
    OnOffFlag.destruct();
    start.destruct();
    interval.destruct();
    TimeStepComp.destruct();
    HarmoComp.destruct();
    HarmoCompByRegions.destruct();
    ResonanceExploration.destruct();
    ResonanceFrequencies.destruct();
    }

private :
 // no private members 
};


//##################################
//
// class of tracers Key Word ....
//
//##################################
class tracers2D_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C start;
  KEY_C interval;
  KEY_C InputFile,InitialisationFile;
  KEY_C TimeIntegration,HorizontalDiffusionMode;
  KEY_C HorizontalDiffusionCoefficient;
  KEY_C SmagorinskyCoefficient;
  KEY_C UpwindCoefficient;
  bool  active;

  tracers2D_C();
  int destruct () {
    OnOffFlag.destruct();
    start.destruct();
    interval.destruct();
    InputFile.destruct();
    InitialisationFile.destruct();
    TimeIntegration.destruct();
    HorizontalDiffusionMode.destruct();
    HorizontalDiffusionCoefficient.destruct();
    SmagorinskyCoefficient.destruct();
    UpwindCoefficient.destruct();
    }

private :
 // no private members 
};

//##################################
//
// class of drifters2D Key Word ....
//
//##################################
class drifters2D_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C start;
  KEY_C interval;
  KEY_C InputFile;
  bool  active;

  drifters2D_C();

  int destruct () {
    OnOffFlag.destruct();
    start.destruct();
    interval.destruct();
    InputFile.destruct();
    }

private :
 // no private members 
};


//##################################
//
// class of dynamics2D Key Word ....
//
//##################################
class dynamics2D_C {
public :
  //constructor puts default value of  variables
  KEY_C WaveEquationForm;
  KEY_C Continuity2DIntegration;
  KEY_C Momentum2DIntegration;
  KEY_C Momentum2DMode;
  KEY_C Advection2Dmode;
  KEY_C AdvectionDiffusionDiscretisation;
  KEY_C SLVDiscretisation;
  KEY_C V2DDiscretisation;
  KEY_C G2DDiscretisation;
  KEY_C WaveScalarProduct;
  KEY_C tau0, theta;
  KEY_C BarotropicW;
  KEY_C HorizontalViscMode,UpwindCoefficient;
  KEY_C AsselinFilteringH, AsselinFilteringU;
  KEY_C AsselinCoefficientC1;
//  KEY_C LateralDiffusion;

  dynamics2D_C();

  int destruct () {
    WaveEquationForm.destruct();
    Continuity2DIntegration.destruct();
    Momentum2DIntegration.destruct();
    Momentum2DMode.destruct();
    Advection2Dmode.destruct();
    AdvectionDiffusionDiscretisation.destruct();
    SLVDiscretisation.destruct();
    V2DDiscretisation.destruct();
    G2DDiscretisation.destruct();
    WaveScalarProduct.destruct();
    tau0.destruct(); theta.destruct();
    BarotropicW.destruct();
    HorizontalViscMode.destruct();UpwindCoefficient.destruct();
    AsselinFilteringH.destruct();AsselinFilteringU.destruct();
    AsselinCoefficientC1.destruct();
    }

private :
 // no private members 
};


//##################################
//
// class of dynamics3D Key Word ....
//
//##################################
class dynamics3D_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C TimeStep;
  KEY_C MomentumIntegration,TracersIntegration;
  KEY_C RhoMode;
  KEY_C TracersMode;
  KEY_C PressureMode;
  KEY_C WDerivation;
  KEY_C TracersDiscretisation,VDiscretisation,WDiscretisation;
  KEY_C TracersAdvectionCentering;
  KEY_C TracerInitialisation; 
  KEY_C InitialisationDirectory;
  KEY_C TFile,SFile,RFile,PFile;
  KEY_C NLayersMax, NUpperHLayers, UpperHLayerLimit;
  KEY_C LevelDistribution;
  KEY_C LevelDisplacement;
  KEY_C TurbulenceClosure;
  KEY_C LateralDiffusion,VerticalDiffusion,TracersVerticalDiffusion;

  dynamics3D_C();
  int destruct () {
   OnOffFlag.destruct();
   TimeStep.destruct();
   MomentumIntegration.destruct();
   TracersIntegration.destruct();
   RhoMode.destruct();
   TracersMode.destruct();
   PressureMode.destruct();
   WDerivation.destruct();
   TracersDiscretisation.destruct();
   VDiscretisation.destruct();
   WDiscretisation.destruct();
   TracersAdvectionCentering.destruct();
   TracerInitialisation.destruct(); 
   InitialisationDirectory.destruct();
   TFile.destruct();
   SFile.destruct();
   RFile.destruct();
   PFile.destruct();
   NLayersMax.destruct();
   NUpperHLayers.destruct();
   UpperHLayerLimit.destruct();
   LevelDistribution.destruct();
   LevelDisplacement.destruct();
   TurbulenceClosure.destruct();
   LateralDiffusion.destruct();VerticalDiffusion.destruct();TracersVerticalDiffusion.destruct();
    }

private :
 // no private members 
};


//##################################
//
// class of tracers3D Key Word ....
//
//##################################
class tracers3D_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C TracersIntegration;
  KEY_C TracersMode;
  KEY_C TracersDiscretisation;
  KEY_C TracersAdvectionCentering;
  KEY_C TracerInitialisation; 
  KEY_C InitialisationDirectory;
  KEY_C TFile,SFile,RFile,PFile;
  KEY_C TVariable;
  KEY_C SVariable;
  KEY_C RVariable;
  KEY_C TurbulenceClosure;
  KEY_C LateralDiffusion,HorizontalDiffusion,VerticalDiffusion;
  KEY_C HorizontalSmooth;
  KEY_C ComputeDiagnostics;

  tracers3D_C();
  int destruct () {
   OnOffFlag.destruct();
   TracersIntegration.destruct();
   TracersMode.destruct();
   TracersDiscretisation.destruct();
   TracersAdvectionCentering.destruct();
   TracerInitialisation.destruct(); 
   InitialisationDirectory.destruct();
   TFile.destruct();SFile.destruct();RFile.destruct();PFile.destruct();
   TVariable.destruct();
   SVariable.destruct();
   RVariable.destruct();
   TurbulenceClosure.destruct();
   LateralDiffusion.destruct();HorizontalDiffusion.destruct();VerticalDiffusion.destruct();
   HorizontalSmooth.destruct();
   ComputeDiagnostics.destruct();
   }

private :
 // no private members 
};

//##################################
//
// class of stream functions Key Word ....
//
//##################################
class streamf_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C ComputeInterval;
  bool  active;

  streamf_C();

  int destruct () {
   OnOffFlag.destruct();
   ComputeInterval.destruct();
    }

private :
 // no private members 
};


//##################################
//
// class of constraint Key Word ....
//
//##################################
class constraint_C {
public :
  //constructor puts default value of  variables
  KEY_C OnOffFlag;
  KEY_C ElevationConstraintFile;
  bool  active;

  constraint_C();

  int destruct () {
    OnOffFlag.destruct();
    ElevationConstraintFile.destruct();
    }

private :
 // no private members 
};


//##################################
//
// class of the tugo configuration ....
//
//##################################

// This is the general class of the tugo configuration
// these class members were initialised by there own constructor ...

class tugo_cfg_C {
public :
  model_C       *model;
  physics_C     *physics;
  dissipation_C *dissipation, *dissipation3D;
  tides_C       *tides;
  atmosphere_C  *atmosphere;
  OceanWaves_C  *OceanWaves;
  boundaries_C  *boundaries;
  ice_C         *ice;
  topography_C     *topography;
  misc_C        *misc;
  fearchive_C   *fearchive;
  structured_C  *structured;
  sample_C      *sample;
  analysis_C    *analysis;
  energy_C      *energy;
  tracers2D_C   *tracers2D;
  drifters2D_C  *drifters2D;
  dynamics2D_C  *dynamics2D;
  dynamics3D_C  *dynamics3D;
  tracers3D_C   *tracers3D;
  streamf_C     *streamf;
  constraint_C  *constraint;

  void initialize() {

    model=0;
    physics=0;
    dissipation=0;
    dissipation3D=0;
    tides=0;
    atmosphere=0;
    OceanWaves=0;
    boundaries=0;
    ice=0;
    topography=0;
    misc=0;
    fearchive=0;
    structured=0;
    sample=0;
    analysis=0;
    energy=0;
    tracers2D=0;
    drifters2D=0;
    dynamics2D=0;
    dynamics3D=0;
    tracers3D=0;
    streamf=0;
    constraint=0;
    };

  tugo_cfg_C() {

    model	 = new model_C;
    physics	 = new physics_C;
    dissipation  = new dissipation_C;
    dissipation3D  = new dissipation_C;
    tides	 = new tides_C;
    atmosphere   = new atmosphere_C;
    OceanWaves   = new OceanWaves_C;
    boundaries   = new boundaries_C;
    ice 	 = new ice_C;
    topography	 = new topography_C;
    misc	 = new misc_C;
    fearchive	 = new fearchive_C;
    structured   = new structured_C;
    sample	 = new sample_C;
    analysis	 = new analysis_C;
    energy	 = new energy_C;
    drifters2D   = new drifters2D_C;
    tracers2D	 = new tracers2D_C;
    tracers3D	 = new tracers3D_C;
    dynamics2D   = new dynamics2D_C;
    dynamics3D   = new dynamics3D_C;
    streamf	 = new streamf_C;
    constraint	 = new constraint_C;
    };


  int destruct () {
    model->destruct();
    physics->destruct();
    dissipation->destruct();
    dissipation3D->destruct();
    tides->destruct();
    atmosphere->destruct();
    OceanWaves->destruct();
    boundaries->destruct();
    ice->destruct();
    topography->destruct();
    misc->destruct();
    fearchive->destruct();
    structured->destruct();
    sample->destruct();
    analysis->destruct();
    energy->destruct();	
    drifters2D->destruct();
    tracers2D->destruct();
    tracers3D->destruct();
    dynamics2D->destruct();
    dynamics3D->destruct();
    streamf->destruct();
    constraint->destruct();

    delete model;
    delete dissipation;
    delete dissipation3D;
    delete tides;
    delete atmosphere;
    delete OceanWaves;
    delete boundaries;
    delete ice;
    delete topography;
    delete misc;
    delete fearchive;
    delete structured;
    delete sample;
    delete analysis;
    delete energy;
    delete drifters2D;
    delete tracers2D;
    delete tracers3D;
    delete dynamics2D;
    delete dynamics3D;
    delete streamf;
    delete constraint;
    };

private :

};


extern void save_partial_config(tugo_cfg_C *tugo_cfg,char *filename);
extern void save_complete_config(tugo_cfg_C *tugo_cfg,char *filename);
extern void hugo_readinput( tugo_cfg_C *tugo_cfg,char *RunControlFile );
extern void import_nl_dot_inp_file(tugo_cfg_C *tugo_cfg,char *filename);
extern int init_config_checks(tugo_cfg_C *tugo_cfg,ostringstream *message,ostringstream *warning,ostringstream *alarm);


#endif

