
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

//  read_init_config.cpp
//
//
//  Toutes ces functions sont utiles pour lire les fichiers de configuration
//  du modele TUGO
//
//  La seule functions qui doit etre appelee depuis l'exterieur
//  est : read_input(tugo_cfg,filemane)
//
//  Thierry LETELLIER LEGOS Toulouse 08Jun2006
//

#include "init_config.hpp" 
#include "keywords.hpp"
#include <string>
#include "exceptions.hpp"
#include <poc-assertions.h>
#include <poc-array.hpp>

extern  bool read_input_model_keyword      ( model_C *,ifstream *in );
extern  bool read_input_physics_keyword    ( physics_C *,ifstream *in );
extern  bool read_input_dissipation_keyword( dissipation_C *,ifstream *in );
extern  bool read_input_tide_keyword       ( tides_C *,ifstream *in );
extern  bool read_input_atmos_keyword      ( atmosphere_C *,ifstream *in );
extern  bool read_input_OceanWaves_keyword ( OceanWaves_C *,ifstream *in );
extern  bool read_input_boundaries_keyword ( boundaries_C *,ifstream *in );
extern  bool read_input_ice_keyword        ( ice_C *,ifstream *in );
extern  bool read_input_topography_keyword    ( topography_C *,ifstream *in );
extern  bool read_input_misc_keyword       ( misc_C *,ifstream *in );
extern  bool read_input_fearch_keyword     ( fearchive_C *,ifstream *in );
extern  bool read_input_structured_keyword ( structured_C *,ifstream *in );
extern  bool read_input_sample_keyword     ( sample_C *,ifstream *in );
extern  bool read_input_analysis_keyword   ( analysis_C *,ifstream *in );
extern  bool read_input_energy_keyword     ( energy_C *,ifstream *in );
extern  bool read_input_drifters2D_keyword ( drifters2D_C *,ifstream *in );
extern  bool read_input_tracers3D_keyword  ( tracers3D_C *,ifstream *in );
extern  bool read_input_tracers_keyword    ( tracers2D_C *,ifstream *in );
extern  bool read_input_dynamics2D_keyword ( dynamics2D_C *,ifstream *in );
extern  bool read_input_dynamics3D_keyword ( dynamics3D_C *,ifstream *in );
extern  bool read_input_streamf_keyword    ( streamf_C *,ifstream *in );
extern  bool read_input_constraint_keyword ( constraint_C *,ifstream *in );


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void exchange_key(KEY_C obsolete, KEY_C *key)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(obsolete.initialised) {
    key->actual_value=obsolete.actual_value;
    obsolete.actual_value="UNDEFINED";
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void solve_deprecated_v01( tugo_cfg_C *tugo_cfg)

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
using namespace Keywords ;
int we_mode,me_mode;
int we_type,me_type;
int u2D_discretisation,z2D_discretisation;
int sproduct2D;

  we_mode = -1;
  if(tugo_cfg->model->WaveEqForm.initialised) {
    if (tugo_cfg->model->WaveEqForm.actual_value      == "P1-P1(lumped)")
      we_mode = 0 ;
    else if (tugo_cfg->model->WaveEqForm.actual_value == "P1-P1(exact)")
      we_mode = 1 ;
    else if (tugo_cfg->model->WaveEqForm.actual_value == "P1-P1")
      we_mode = 1 ;
    else if (tugo_cfg->model->WaveEqForm.actual_value == "P1-NCP1")
      we_mode = 2 ;
    else if (tugo_cfg->model->WaveEqForm.actual_value == "WE-explicit")
      we_mode = 3 ;
    else if (tugo_cfg->model->WaveEqForm.actual_value == "WE-implicit")
      we_mode = 4 ;
    else 
      check_error(-1, "input error", __LINE__, __FILE__, 1);
    }

  tugo_cfg->model->WaveEqForm.actual_value="UNDEFINED";

  me_mode = -1;
  if(tugo_cfg->model->TwoDMometumEqForm.initialised) {
    if (tugo_cfg->model->TwoDMometumEqForm.actual_value      == "P1-P1")
      me_mode = 0 ;
    else if (tugo_cfg->model->TwoDMometumEqForm.actual_value == "P1-NCP1-fake")
      me_mode = 1 ;
    else if (tugo_cfg->model->TwoDMometumEqForm.actual_value == "P1-NCP1")
      me_mode = 2 ;
    else 
      check_error(-1, "input error", __LINE__, __FILE__, 1);
    }

  tugo_cfg->model->TwoDMometumEqForm.actual_value="UNDEFINED";

  switch (we_mode) {
    case 0:
      tugo_cfg->dynamics2D->WaveEquationForm.actual_value   = KEY_GWE_CLX;
      tugo_cfg->dynamics2D->Continuity2DIntegration.actual_value   = KEY_UNSPECIFIED;
      tugo_cfg->dynamics2D->SLVDiscretisation.actual_value  = KEY_LGP1;
      tugo_cfg->dynamics2D->V2DDiscretisation.actual_value  = KEY_LGP1;
      tugo_cfg->dynamics2D->WaveScalarProduct.actual_value  = KEY_NQUAD;
      break;

    case 1:
      tugo_cfg->dynamics2D->WaveEquationForm.actual_value   = KEY_GWE_CLX;
      tugo_cfg->dynamics2D->Continuity2DIntegration.actual_value   = KEY_UNSPECIFIED;
      tugo_cfg->dynamics2D->SLVDiscretisation.actual_value  = KEY_LGP1_CG;
      tugo_cfg->dynamics2D->V2DDiscretisation.actual_value  = KEY_LGP1_CG;
      tugo_cfg->dynamics2D->WaveScalarProduct.actual_value  = KEY_INTGL;
      break;

    case 2:
      tugo_cfg->dynamics2D->WaveEquationForm.actual_value   = KEY_GWE_CLX;
      tugo_cfg->dynamics2D->Continuity2DIntegration.actual_value   = KEY_UNSPECIFIED;
      tugo_cfg->dynamics2D->SLVDiscretisation.actual_value  = KEY_LGP1_CG;
      tugo_cfg->dynamics2D->V2DDiscretisation.actual_value  = KEY_NCP1_CG;
      tugo_cfg->dynamics2D->WaveScalarProduct.actual_value  = KEY_INTGL;
      break;

    case 3:
      tugo_cfg->dynamics2D->WaveEquationForm.actual_value   = KEY_WE_EXPLICIT;
      tugo_cfg->dynamics2D->Continuity2DIntegration.actual_value   = KEY_LEAPFROG;
      tugo_cfg->dynamics2D->SLVDiscretisation.actual_value  = KEY_LGP1_CG;
      tugo_cfg->dynamics2D->V2DDiscretisation.actual_value  = KEY_NCP1_CG;
      tugo_cfg->dynamics2D->WaveScalarProduct.actual_value  = KEY_INTGL;
      break;

    case 4:
      tugo_cfg->dynamics2D->WaveEquationForm.actual_value   = KEY_WE_IMPLICIT;
      tugo_cfg->dynamics2D->Continuity2DIntegration.actual_value   = KEY_LEAPFROG;
      tugo_cfg->dynamics2D->SLVDiscretisation.actual_value  = KEY_LGP1_CG;
      tugo_cfg->dynamics2D->V2DDiscretisation.actual_value  = KEY_NCP1_CG;
      tugo_cfg->dynamics2D->WaveScalarProduct.actual_value  = KEY_INTGL;
      break;

//     default:
//       check_error(-1, "illegal  mode error", __LINE__, __FILE__, 1);
//       break;
    }

  switch (me_mode) {
    case 0:
      tugo_cfg->dynamics2D->V2DDiscretisation.actual_value  = KEY_LGP1;
      break;

    case 1:
      tugo_cfg->dynamics2D->V2DDiscretisation.actual_value  = KEY_P1NCP1FAKE;
      break;

    case 2:
      tugo_cfg->dynamics2D->V2DDiscretisation.actual_value  = KEY_NCP1;
      break;
    }

  if(tugo_cfg->model->theta.initialised) {
    tugo_cfg->dynamics2D->theta.actual_value=tugo_cfg->model->theta.actual_value;
    tugo_cfg->model->theta.actual_value="UNDEFINED";
    }

  if(tugo_cfg->model->tau0.initialised) {
    tugo_cfg->dynamics2D->tau0.actual_value=tugo_cfg->model->tau0.actual_value;
    tugo_cfg->model->tau0.actual_value="UNDEFINED";
    }

  if(tugo_cfg->misc->BarotropicW.initialised) {
    tugo_cfg->dynamics2D->BarotropicW.actual_value=tugo_cfg->misc->BarotropicW.actual_value;
    tugo_cfg->misc->BarotropicW.actual_value="UNDEFINED";
    }



}
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void solve_deprecated_v02( tugo_cfg_C *tugo_cfg)

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
using namespace Keywords ;

  exchange_key(tugo_cfg->dynamics3D->TracersIntegration,&tugo_cfg->tracers3D->TracersIntegration);
  exchange_key(tugo_cfg->dynamics3D->TracersMode,&tugo_cfg->tracers3D->TracersMode);
  exchange_key(tugo_cfg->dynamics3D->TracersDiscretisation,&tugo_cfg->tracers3D->TracersDiscretisation); 
  exchange_key(tugo_cfg->dynamics3D->TracersAdvectionCentering,&tugo_cfg->tracers3D->TracersAdvectionCentering);
  exchange_key(tugo_cfg->dynamics3D->TracerInitialisation,&tugo_cfg->tracers3D->TracerInitialisation);
  exchange_key(tugo_cfg->dynamics3D->InitialisationDirectory,&tugo_cfg->tracers3D->InitialisationDirectory);
  exchange_key(tugo_cfg->dynamics3D->TFile,&tugo_cfg->tracers3D->TFile);
  exchange_key(tugo_cfg->dynamics3D->SFile,&tugo_cfg->tracers3D->SFile);
  exchange_key(tugo_cfg->dynamics3D->RFile,&tugo_cfg->tracers3D->RFile);
  exchange_key(tugo_cfg->dynamics3D->PFile,&tugo_cfg->tracers3D->PFile);
  exchange_key(tugo_cfg->dynamics3D->TurbulenceClosure,&tugo_cfg->tracers3D->TurbulenceClosure);
  exchange_key(tugo_cfg->dynamics3D->TracersVerticalDiffusion,&tugo_cfg->tracers3D->VerticalDiffusion);

  exchange_key(tugo_cfg->structured->ArchiveInterval,&tugo_cfg->structured->Interval);
  exchange_key(tugo_cfg->structured->ArchiveStart,&tugo_cfg->structured->Start);

  exchange_key(tugo_cfg->fearchive->ArchOverAppendFlag,&tugo_cfg->fearchive->OverwriteAppendFlag);

  exchange_key(tugo_cfg->dissipation->KarmanRugosity,&tugo_cfg->dissipation->RoughnessLength);
  
  exchange_key(tugo_cfg->model->TopoFile,&tugo_cfg->topography->TopoFile);
  exchange_key(tugo_cfg->model->SlopeFile,&tugo_cfg->topography->SlopeFile);
  exchange_key(tugo_cfg->model->ZminFile,&tugo_cfg->topography->ZminFile);
  exchange_key(tugo_cfg->model->RelativeMinDepth,&tugo_cfg->topography->RelativeMinDepth);
  exchange_key(tugo_cfg->model->MinDepth,&tugo_cfg->topography->MinDepth);

//   exchange_key(tugo_cfg->model->,&tugo_cfg->topography->);
// 
//   TopoFile.deprecated=true;
//   SlopeFile.deprecated=true;
//   ZminFile.deprecated=true;
//   RelativeMinDepth.deprecated=true;
//   MinDepth.deprecated=true;
  
}
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void hugo_readinput( tugo_cfg_C *tugo_cfg,char *RunControlFile )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ifstream in;
  string s;
  int pos;
  void *keyword;
  bool the_test;
  int status;

  in.open(RunControlFile, ifstream::in);

  if(in.fail()) {
    check_error(-1, "cannot open T-UGO input file", __LINE__, __FILE__, 1);
    }

  the_test==false;

  while( (getline(in,s)) ) {
    pos=s.find("%");
    if(pos!=-1) { //I found a % in the line then I simply copy the line in the echo file
      if(pos!=0) cout << "find a % at position " << pos << " in this line : \"" << s << "\" ----> presume that this is a comment" << endl << flush;
      continue;
      }
    pos=s.find("#");
    if(pos!=-1) { //I found a keyword then I have to treat the next keys
      if( s.find("model")==pos+1 )          {the_test=read_input_model_keyword(tugo_cfg->model , &in);continue; }
      if( s.find("physics")==pos+1 )        {the_test=read_input_physics_keyword(tugo_cfg->physics , &in);continue; }
      if( s.find("dissipation3D")==pos+1 )  { 
        the_test=read_input_dissipation_keyword(tugo_cfg->dissipation3D , &in);continue;
        }
      if( s.find("dissipation")==pos+1 )    {the_test=read_input_dissipation_keyword(tugo_cfg->dissipation , &in);continue; }
      if( s.find("tide")==pos+1 )           {the_test=read_input_tide_keyword(tugo_cfg->tides , &in);continue; }
      if( s.find("atmosphere")==pos+1 )     {the_test=read_input_atmos_keyword(tugo_cfg->atmosphere , &in);continue; }
      if( s.find("OceanWaves")==pos+1 )     {
        the_test=read_input_OceanWaves_keyword(tugo_cfg->OceanWaves , &in);continue;
        }
      if( s.find("boundaries")==pos+1 )     {the_test=read_input_boundaries_keyword(tugo_cfg->boundaries , &in);continue; }
      if( s.find("ice")==pos+1 )            {the_test=read_input_ice_keyword(tugo_cfg->ice ,               &in);continue; }
      if( s.find("miscellaneous")==pos+1 )  {the_test=read_input_misc_keyword(tugo_cfg->misc ,             &in);continue; }
      if( s.find("topography")==pos+1 )     {the_test=read_input_topography_keyword(tugo_cfg->topography,  &in);continue; }
/**----------------------------------------------------------------------------
      redundant, to stay compatible with earlier format*/
      if( s.find("tsunami")==pos+1 )        {
	the_test=read_input_topography_keyword(tugo_cfg->topography,  &in);continue;
        }
      if( s.find("fe_archive")==pos+1 )     {
        the_test=read_input_fearch_keyword(tugo_cfg->fearchive , &in);
        continue;
        }
      if( s.find("structured_archive")==pos+1 ) {
//        tugo_cfg->structured->Start=tugo_cfg->fearchive->ArchStart;
//        tugo_cfg->structured->Interval=tugo_cfg->fearchive->ArchInterval;
        the_test=read_input_structured_keyword(tugo_cfg->structured , &in);
        continue;
        }
      if( s.find("sample_points")==pos+1 )   {the_test=read_input_sample_keyword(tugo_cfg->sample ,         &in);continue; }
      if( s.find("analysis")==pos+1 )        {the_test=read_input_analysis_keyword(tugo_cfg->analysis ,     &in);continue; }
      if( s.find("energy_budget")==pos+1 )   {the_test=read_input_energy_keyword(tugo_cfg->energy ,         &in);continue; }
      if( s.find("drifters2D")==pos+1 )      {the_test=read_input_drifters2D_keyword(tugo_cfg->drifters2D , &in);continue; }
/**----------------------------------------------------------------------------
      redundant, to stay compatible with earlier format*/
      if( s.find("drifters")==pos+1 )        {the_test=read_input_drifters2D_keyword(tugo_cfg->drifters2D , &in);continue; }
      if( s.find("tracers3D")==pos+1 )       {
        the_test=read_input_tracers3D_keyword(tugo_cfg->tracers3D , &in);continue;
        }
      if( s.find("tracers2D")==pos+1 )       {the_test=read_input_tracers_keyword(tugo_cfg->tracers2D , &in);continue; }
/**----------------------------------------------------------------------------
      redundant, to stay compatible with earlier format*/
      if( s.find("tracers")==pos+1 )       {
        the_test=read_input_tracers_keyword(tugo_cfg->tracers2D , &in);continue;
        }
      if( s.find("dynamics2D")==pos+1 )      {the_test=read_input_dynamics2D_keyword(tugo_cfg->dynamics2D , &in);continue; }
      if( s.find("dynamics3D")==pos+1 )      {the_test=read_input_dynamics3D_keyword(tugo_cfg->dynamics3D , &in);continue; }
      if( s.find("streamf")==pos+1 )         {the_test=read_input_streamf_keyword(tugo_cfg->streamf ,       &in);continue; }
      if( s.find("constraints")==pos+1 )     {the_test=read_input_constraint_keyword(tugo_cfg->constraint , &in);continue; }
      }
    cout << __FUNCTION__ <<": can't understand this line \n--> " <<  s <<  endl << flush;
    }

  solve_deprecated_v01(tugo_cfg);
  solve_deprecated_v02(tugo_cfg);

}

//
// the next functions sould not be used alone ....
//

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static  void get_char_value(string s,KEY_C *key)

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string comments;
  bool recommence;

/*----------------------------------------------------------------------
  tout d'abord je cherche le signe =  */
  pos=s.find("=");
 
/*----------------------------------------------------------------------
  et j'enleve tout ce qui est devant*/
  s=s.substr(pos+1);
  
  pos=s.find("<default>");
  if(pos>0) {
    s=s.substr(pos+9);
    }

/*----------------------------------------------------------------------
  recherche des commentaires
  copy dans la chaine comments*/
  comments.assign(s);
  pos=s.find("//");
  comments=comments.substr(pos+2);

/*----------------------------------------------------------------------
  maintenant j'enleve les commentaires dans s */
  //s = s.erase(s.size() - comments.size() - 2, comments.size() + 2) ;
  s = string(s, 0, pos - 1) ;
/*----------------------------------------------------------------------
  maintenant j'enleve les espaces a la fin*/
  recommence=true;
  while(recommence){
    if( s.rfind(" ")==s.size()-1 )
      {
      s=s.erase(s.size()-1,1 );
      recommence=true;
      }
    else recommence =false;
    }

/*----------------------------------------------------------------------
  enfin j'enleve les espaces au debut*/
  recommence=true;
  while(recommence){
    if( s.find(" ")==0 ) {
      s=s.erase(0,1 );
      recommence=true;
      }
    else recommence =false;
    }

/*----------------------------------------------------------------------
  il me reste a remplir la key*/

  //sprintf(key->actual_value,"%s",s.c_str() );
  key->actual_value = s ;
  if(comments!="some comments") {
    sprintf(key->comments,"%s",comments.c_str() );
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static bool identify_key(string s,KEY_C *key)

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string keyword(key->keyword); 

  keyword.append(" ");

/**----------------------------------------------------------------------
  danger if 1 key is entirely part of another*/
  pos=s.find(keyword); 
  if( (pos!=-1)&&(pos<10) ){
    if(key->deprecated) {
      cout << __FUNCTION__ <<": deprecated keyword \n--> " <<  key->keyword <<  endl << flush;
      }
    if(pos>1) {
      if(s[pos-1]!=' ') return(false);
      }
    get_char_value(s,key);
    key->initialised=true;
    return(true);}
  else return(false);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_model_keyword( model_C *model,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

//  cout <<__FUNCTION__<< "\n" << flush;
  
  while(getline(*in,s) ) {
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1){return(true);}

    if(identify_key(s,&model->rootname) )   continue;
    if(identify_key(s,&model->MeshFormat) ) continue;
    if(identify_key(s,&model->MeshType) ) {
/**----------------------------------------------------------------------------
      for backward compatibility reasons */
      if(model->MeshType.actual_value=="nei") {
        printf("obsolete input, patch applied; please re-process later through tugo-gui\n");
        model->MeshType.actual_value.assign("TRIANGULAR");
        model->MeshFormat.actual_value.assign("nei");
        }
      continue;
      }
    if(identify_key(s,&model->MeshMeta) ) continue;
    if(identify_key(s,&model->MeshFile) ) continue;
    if(identify_key(s,&model->BelFile) ) {
      continue;
      }
    if(identify_key(s,&model->TopoFile) )             continue;
    if(identify_key(s,&model->SlopeFile) )            continue;
    if(identify_key(s,&model->ZminFile) )             continue;
    if(identify_key(s,&model->RelativeMinDepth) )     continue;
    if(identify_key(s,&model->MinDepth) )             continue;
    if(identify_key(s,&model->OrigineTime) )          continue;
    if(identify_key(s,&model->SubTimeStep) )          continue;
    if(identify_key(s,&model->TimeStep) )             continue;
    if(identify_key(s,&model->RunDuration) )          continue;
    if(identify_key(s,&model->SpinUpDuration) )       continue;
    if(identify_key(s,&model->SpinUpMode) )           continue;
    if(identify_key(s,&model->RestartFile) )          continue;
    if(identify_key(s,&model->RestartFormat) )        continue;
    if(identify_key(s,&model->ContinuationFile) )     continue;
    if(identify_key(s,&model->ContinuationFormat) )   continue;
    if(identify_key(s,&model->ContinuationInterval) ) continue;
    if(identify_key(s,&model->AutomaticInterval) )    continue;
    if(identify_key(s,&model->OutputPath) )           continue;
    if(identify_key(s,&model->LinearSolver) )         continue;
    if(identify_key(s,&model->SubLinearSolver) )      continue;
    if(identify_key(s,&model->SolverMode) )           continue;
    if(identify_key(s,&model->mode) )                 continue;
    if(identify_key(s,&model->compressibility) )      continue;
    if(identify_key(s,&model->WaveEqForm) )           continue;
    if(identify_key(s,&model->TwoDMometumEqForm) )    continue;
    if(identify_key(s,&model->tau0) )                 continue;
    if(identify_key(s,&model->theta) )                continue;
    if(identify_key(s,&model->OPENMP_nCPUsMax) )      continue;
    if(identify_key(s,&model->MetisPath) )            continue;
    if(identify_key(s,&model->ThreeDPressure) )       continue;
    if(identify_key(s,&model->EchoFile) )             continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_physics_keyword( physics_C *physics,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

//  cout <<__FUNCTION__<< "\n" << flush;
  
  while(getline(*in,s) ) {
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1){return(true);}

    if(identify_key(s,&physics->gravity_mode) )        continue;
    if(identify_key(s,&physics->anelastic_body_tide) ) continue;
    if(identify_key(s,&physics->gravity_adjustment) )  continue;
    if(identify_key(s,&physics->coriolis_mode) )       continue;
    if(identify_key(s,&physics->coriolis_latitude) )   continue;
    if(identify_key(s,&physics->water_viscosity) )     continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_dynamics2D_keyword( dynamics2D_C *dynamics2D,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    //energy->compute=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1){return(true);}

    if(identify_key(s,&dynamics2D->WaveEquationForm) )                 continue;
    if(identify_key(s,&dynamics2D->Continuity2DIntegration) )          continue;
    if(identify_key(s,&dynamics2D->Momentum2DIntegration) )            continue;
    if(identify_key(s,&dynamics2D->Momentum2DMode) )                   continue;
    if(identify_key(s,&dynamics2D->Advection2Dmode) )                  continue;
    if(identify_key(s,&dynamics2D->AdvectionDiffusionDiscretisation) ) continue;
    if(identify_key(s,&dynamics2D->SLVDiscretisation) )                continue;
    if(identify_key(s,&dynamics2D->V2DDiscretisation) )                continue;
    if(identify_key(s,&dynamics2D->G2DDiscretisation) )                continue;
    if(identify_key(s,&dynamics2D->WaveScalarProduct) )                continue;
    if(identify_key(s,&dynamics2D->tau0) )                             continue;
    if(identify_key(s,&dynamics2D->theta) )                            continue;
    if(identify_key(s,&dynamics2D->BarotropicW) )                      continue;
    if(identify_key(s,&dynamics2D->HorizontalViscMode) )               continue;
    if(identify_key(s,&dynamics2D->UpwindCoefficient) )                continue;
    if(identify_key(s,&dynamics2D->AsselinFilteringH) )                continue;
    if(identify_key(s,&dynamics2D->AsselinFilteringU) )                continue;
    if(identify_key(s,&dynamics2D->AsselinCoefficientC1) )                continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_dynamics3D_keyword( dynamics3D_C *dynamics3D,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    //energy->compute=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1){return(true);}

    if(identify_key(s,&dynamics3D->OnOffFlag) ) continue;
    if(identify_key(s,&dynamics3D->TimeStep) )  continue;
    if(identify_key(s,&dynamics3D->MomentumIntegration) ) continue;
    if(identify_key(s,&dynamics3D->TracersIntegration) )  continue;
    if(identify_key(s,&dynamics3D->RhoMode) )      continue;
    if(identify_key(s,&dynamics3D->TracersMode) )  continue;
    if(identify_key(s,&dynamics3D->PressureMode) ) continue;
    if(identify_key(s,&dynamics3D->WDerivation) )  continue;
    if(identify_key(s,&dynamics3D->VDiscretisation) ) continue;
    if(identify_key(s,&dynamics3D->WDiscretisation) ) continue;
    
    if(identify_key(s,&dynamics3D->TracersDiscretisation) )     continue;
    if(identify_key(s,&dynamics3D->TracersAdvectionCentering) ) continue;
    if(identify_key(s,&dynamics3D->TracerInitialisation) )      continue;
    if(identify_key(s,&dynamics3D->InitialisationDirectory) )   continue;
    
    if(identify_key(s,&dynamics3D->TFile) ) continue;
    if(identify_key(s,&dynamics3D->SFile) ) continue;
    if(identify_key(s,&dynamics3D->RFile) ) continue;
    if(identify_key(s,&dynamics3D->PFile) ) continue;
    
    if(identify_key(s,&dynamics3D->NLayersMax) )        continue;
    if(identify_key(s,&dynamics3D->NUpperHLayers) )     continue;
    if(identify_key(s,&dynamics3D->UpperHLayerLimit) )  continue;
    if(identify_key(s,&dynamics3D->LevelDistribution) ) continue;
    
    if(identify_key(s,&dynamics3D->LevelDisplacement) ) continue;
    if(identify_key(s,&dynamics3D->TurbulenceClosure) ) continue;
    if(identify_key(s,&dynamics3D->LateralDiffusion) )  continue;
//    if(identify_key(s,&dynamics3D->HorizontalDiffusion) )  continue;
    if(identify_key(s,&dynamics3D->VerticalDiffusion) ) continue;
    if(identify_key(s,&dynamics3D->TracersVerticalDiffusion) ) continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_dissipation_keyword( dissipation_C *dissipation,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1){return(true);}

    if(identify_key(s,&dissipation->MinBackGroundSpd) )       continue;
    if(identify_key(s,&dissipation->BKGPolygons) )            continue;
    if(identify_key(s,&dissipation->BKGValues) )              continue;
    if(identify_key(s,&dissipation->BackgroundVelocityFile) ) continue;
    if(identify_key(s,&dissipation->IceShelvesFile) )         continue;
    if(identify_key(s,&dissipation->IceShelvesCoefficient) )  continue;
    if(identify_key(s,&dissipation->BottomFrictionType) )     continue;
    if(identify_key(s,&dissipation->BottomRugosityFile) )     continue;
    if(identify_key(s,&dissipation->LinearFrictionCoeff) )    continue;
    if(identify_key(s,&dissipation->QuadraticFrictionCoeff) ) continue;
    if(identify_key(s,&dissipation->QFCPolygons) )            continue;
    if(identify_key(s,&dissipation->QFCValues) )              continue;
    if(identify_key(s,&dissipation->KarmanRugosity) )         continue;
    if(identify_key(s,&dissipation->RoughnessLength) )        continue;
    if(identify_key(s,&dissipation->MinQuadraticFrictionCoeff) ) continue;
    if(identify_key(s,&dissipation->QFCMPolygons) )           continue;
    if(identify_key(s,&dissipation->QFCMValues) )             continue;
    if(identify_key(s,&dissipation->BottomFrictionCoeffFile) ) continue;
    if(identify_key(s,&dissipation->RugosityPolygons) )       continue;
    if(identify_key(s,&dissipation->RugosityValues) )         continue;
    if(identify_key(s,&dissipation->RugosityValuesbyRegions) ) continue;
    if(identify_key(s,&dissipation->InternDragAlgorithm) )    continue;
    if(identify_key(s,&dissipation->InternDragSlope) )        continue;
    if(identify_key(s,&dissipation->InternDragHmin) )         continue;
    if(identify_key(s,&dissipation->InternDragHmax) )         continue;
    if(identify_key(s,&dissipation->IWD01Polygons) )          continue;
    if(identify_key(s,&dissipation->IWD01Values) )            continue;
    if(identify_key(s,&dissipation->IWD01ValuesByRegions) )   continue;
    if(identify_key(s,&dissipation->InternDragRugosity) )     continue;
    if(identify_key(s,&dissipation->MixedLayerCoeff) )        continue;
    if(identify_key(s,&dissipation->MinhorizontalVisc) )      continue;
    if(identify_key(s,&dissipation->SmagorinskyCoefficient) ) continue;
    if(identify_key(s,&dissipation->HorizontalViscMode) )     continue;
    if(identify_key(s,&dissipation->BruntVassalaValue) )      continue;
    if(identify_key(s,&dissipation->BruntVassalaFile) )       continue;
    if(identify_key(s,&dissipation->ShearDragFlag) )          continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout << "can not understand this line" << s  << "\n" << flush;continue;
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_tide_keyword( tides_C *tide,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    //cout << s <<flush;
    tide->active=true;

    pos=s.find("##"); if(pos!=-1) goto finished;

    if(identify_key(s,&tide->OnOffFlag) )               continue;
    if(identify_key(s,&tide->AstronomicPotentialFlag) ) continue;
    if(identify_key(s,&tide->LSAFlag) )                 continue;
    if(identify_key(s,&tide->AtlasDir) )                continue;
    if(identify_key(s,&tide->AtlasConvention) )         continue;
    if(identify_key(s,&tide->LSADir) )                  continue;
    if(identify_key(s,&tide->LSAConvention) )           continue;
    if(identify_key(s,&tide->PressureFlag) )            continue;
    if(identify_key(s,&tide->PressureDir) )             continue;
    if(identify_key(s,&tide->PressureConvention) )      continue;
    if(identify_key(s,&tide->admittance) )              continue;
    if(identify_key(s,&tide->equilibrium) )             continue;
    if(identify_key(s,&tide->NodalCorrection) )         continue;
    if(identify_key(s,&tide->LSACoeff) )                continue;
    if(identify_key(s,&tide->DeformationLoveNum) )      continue;
    if(identify_key(s,&tide->BoundTideFlag) )           continue;
    if(identify_key(s,&tide->BoundTideFile) )           continue;
    if(identify_key(s,&tide->Spectral2DRun) )           continue;
    if(identify_key(s,&tide->SpectralSolver) )          continue;
    if(identify_key(s,&tide->DominantWave_1) )          continue;
    if(identify_key(s,&tide->DominantWave_2) )          continue;
    if(identify_key(s,&tide->PriorSolution) )          continue;
    if(identify_key(s,&tide->SpectralPaire) )           continue;
    if(identify_key(s,&tide->Spectral2DMaxIteration) )  continue;
    if(identify_key(s,&tide->CompoundMaxIteration) )    continue;
    if(identify_key(s,&tide->SpecificData) )            continue;
    if(identify_key(s,&tide->TopoFile) )                continue;
    if(identify_key(s,&tide->SlopeFile) )               continue;
    if(identify_key(s,&tide->SpectralRecycle) )         continue;
    if(identify_key(s,&tide->ReducedArchive) )          continue;
    if(identify_key(s,&tide->SpectralArchiveLevel) )    continue;
    if(identify_key(s,&tide->Spectral3DRun) )           continue;
    if(identify_key(s,&tide->Spectral3DOBCs) )          continue;
    if(identify_key(s,&tide->VDiffusionDuration) )      continue;
    if(identify_key(s,&tide->VDiffusionTimeStep) )      continue;
    if(identify_key(s,&tide->Spectral3DFirstIteration) )  continue;
    if(identify_key(s,&tide->Spectral3DMaxIteration) )  continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;
    }

finished:
  if(tide->OnOffFlag.actual_value == "FALSE") {
    tide->active=false;
    }
  return(true);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_atmos_keyword( atmosphere_C *atmosphere,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) ) {
    //cout << s <<flush;
    atmosphere->active=true;
    pos=s.find("##"); if(pos!=-1) goto finished;

    if(identify_key(s,&atmosphere->OnOffFlag) )             continue;
    if(identify_key(s,&atmosphere->PressureFlag) )          continue;
    if(identify_key(s,&atmosphere->PressureDir) ) {
//      printf("%s\n",s.c_str());
      continue;
      }
    if(identify_key(s,&atmosphere->PressureConvention) )    continue;
    if(identify_key(s,&atmosphere->WindFlag) )              continue;
    if(identify_key(s,&atmosphere->LongWave) )              continue;
    if(identify_key(s,&atmosphere->ShortWave) )             continue;
    if(identify_key(s,&atmosphere->AtmosDirectory) )        continue;
    if(identify_key(s,&atmosphere->FileFormat) )            continue;
    if(identify_key(s,&atmosphere->FileEndian) )            continue;
    if(identify_key(s,&atmosphere->FileConvention) )        continue;
    if(identify_key(s,&atmosphere->MaskFile) )              continue;
    if(identify_key(s,&atmosphere->ApplyMask) )             continue;
    if(identify_key(s,&atmosphere->KeepMeanPressureFlag) )  continue;
    if(identify_key(s,&atmosphere->KeepMeanWindFlag) )      continue;
    if(identify_key(s,&atmosphere->TimeInterp) )            continue;
    if(identify_key(s,&atmosphere->BulkFormula) )           continue;
    if(identify_key(s,&atmosphere->WindStress) )           continue;
    if(identify_key(s,&atmosphere->BoundInversebarometer) ) continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
finished:
  atmosphere->active=(atmosphere->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_OceanWaves_keyword(OceanWaves_C *keystruc,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) ) {
    //cout << s <<flush;
    keystruc->active=true;
    pos=s.find("##"); if(pos!=-1) goto finished;

    if(identify_key(s,&keystruc->OnOffFlag) )         continue;
    if(identify_key(s,&keystruc->Parameterisation) )  continue;
    if(identify_key(s,&keystruc->DragFlag) )          continue;
    if(identify_key(s,&keystruc->FrictionFlag) )      continue;
    if(identify_key(s,&keystruc->TimeInterp) )        continue;
    if(identify_key(s,&keystruc->Discretisation) )    continue;
    if(identify_key(s,&keystruc->FileDirectory) )     continue;
    if(identify_key(s,&keystruc->FileFormat) )        continue;
    if(identify_key(s,&keystruc->FileConvention) )    continue;
    if(identify_key(s,&keystruc->MaskFile) )          continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
finished:
  keystruc->active=(keystruc->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_boundaries_keyword( boundaries_C *bound,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1){return(true);}    

    if(identify_key(s,&bound->type) )                 continue;
    if(identify_key(s,&bound->mode) )                 continue;
    if(identify_key(s,&bound->RelaxTime) )            continue;
    if(identify_key(s,&bound->ArchivedBCUse) )        continue;
    if(identify_key(s,&bound->ArchivedBCDirectory) )  continue;
    if(identify_key(s,&bound->ArchivedBCFormat) )     continue;
    if(identify_key(s,&bound->ArchivedBCRoot) )       continue;
    if(identify_key(s,&bound->ArchivedBCUnits) )      continue;
    if(identify_key(s,&bound->ArchivedBCTransport) )  continue;
    if(identify_key(s,&bound->PeriodicFile) )         continue;
    if(identify_key(s,&bound->SpongeFlag) )           continue;
    if(identify_key(s,&bound->BufferZoneFile) )       continue;
    if(identify_key(s,&bound->PMLfactor) )            continue;
    if(identify_key(s,&bound->RiversFlag) )           continue;
    if(identify_key(s,&bound->RiversFile) )           continue;
    if(identify_key(s,&bound->RiversTimeTemplate) )   continue;
    if(identify_key(s,&bound->SolidCondition) )       continue;
    if(identify_key(s,&bound->NoSlipSize) )           continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_ice_keyword( ice_C *ice,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    //cout << s <<flush;
    ice->active=true;
    pos=s.find("##"); if(pos!=-1) goto finished;

    if(identify_key(s,&ice->OnOffFlag) )         continue;
    if(identify_key(s,&ice->TimeStep) )          continue;
    if(identify_key(s,&ice->MomentumSolver) )    continue;
    if(identify_key(s,&ice->Rheology) )          continue;
    if(identify_key(s,&ice->IceCoverFile) )      continue;
    if(identify_key(s,&ice->IceElasticityFile) ) continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
finished:
  ice->active=(ice->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_topography_keyword( topography_C *topography,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1)goto finished;  

    if(identify_key(s,&topography->OnOffFlag) )               continue;
    if(identify_key(s,&topography->BottomMotionFlag) )        continue;
    if(identify_key(s,&topography->BottomMotionDirectory) )   continue;
    if(identify_key(s,&topography->BottomMotionFileConvent) ) continue;
    
    if(identify_key(s,&topography->TopoFile) )             continue;
    if(identify_key(s,&topography->SlopeFile) )            continue;
    if(identify_key(s,&topography->ZminFile) )             continue;
    if(identify_key(s,&topography->RelativeMinDepth) )     continue;
    if(identify_key(s,&topography->MinDepth) )             continue;
    if(identify_key(s,&topography->DepthOffsetPolygons) )  continue;
    if(identify_key(s,&topography->DepthOffsetValues) )    continue;
    if(identify_key(s,&topography->MeanLevelFile) )        continue;
    if(identify_key(s,&topography->DepthScaleFactor) )    continue;
    if(identify_key(s,&topography->DryingMinDepth) )    continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
finished:
  topography->BottomMotionActive=(topography->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  topography->active=true;
  return(true);
}
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_misc_keyword( misc_C *misc,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1){return(true);}      

    if(identify_key(s,&misc->NonTidalLoadFlag) )   continue;
    if(identify_key(s,&misc->SmoothingFlag) )      continue;
    if(identify_key(s,&misc->checks) )             continue;
    if(identify_key(s,&misc->Check_Wconsistency) ) continue;
    if(identify_key(s,&misc->CompatibilityScale) ) continue;
    if(identify_key(s,&misc->BarotropicW) )        continue;
    if(identify_key(s,&misc->coupling_2D_3D) )     continue;
    if(identify_key(s,&misc->StabilityControl) )   continue;
    if(identify_key(s,&misc->StabilityCheck) )     continue;
    if(identify_key(s,&misc->MaxDeltaElevation) )  continue;
    if(identify_key(s,&misc->MaxDeltaCurrents) )   continue;
    if(identify_key(s,&misc->MaxDeltaPeriod) )     continue;
    if(identify_key(s,&misc->nSubCycles) )         continue;
    if(identify_key(s,&misc->SubcyclingPersistence) )   continue;
    if(identify_key(s,&misc->SubcyclingTemporization) ) continue;
    if(identify_key(s,&misc->SubcyclingExtension) )     continue;
    if(identify_key(s,&misc->SubcyclingLockLimit) )     continue;
    if(identify_key(s,&misc->SubcyclingCFLratio) )     continue;
    if(identify_key(s,&misc->SmoothZField) )        continue;
    if(identify_key(s,&misc->SmoothUField) )        continue;
    if(identify_key(s,&misc->SmoothAdvection) )     continue;
    if(identify_key(s,&misc->SmoothDiffusion) )     continue;
    if(identify_key(s,&misc->FlushTracking) ) {
      continue;
      };
    if(identify_key(s,&misc->ToponymsFile) )     continue;
    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_fearch_keyword( fearchive_C *fearch,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    fearch->active=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1)goto finished;

    if(identify_key(s,&fearch->OnOffFlag) ) continue;
    if(identify_key(s,&fearch->StateFile) ) continue;
    if(identify_key(s,&fearch->VectorFile) ) continue;
    if(identify_key(s,&fearch->ArchFormat) ) continue;
    if(identify_key(s,&fearch->ArchConvention) ) continue;
    if(identify_key(s,&fearch->ArchStart) ) continue;
    if(identify_key(s,&fearch->ArchInterval) ) continue;
    if(identify_key(s,&fearch->ArchBits) ) continue;
    if(identify_key(s,&fearch->Arch_huv_Scale) ) continue;
    if(identify_key(s,&fearch->Arch_pW_Scale) ) continue;
    if(identify_key(s,&fearch->ArchOverAppendFlag) ) continue;
    if(identify_key(s,&fearch->OverwriteAppendFlag) ) continue;

    if(identify_key(s,&fearch->z_flag) ) continue;
    if(identify_key(s,&fearch->ubar_flag) ) continue;
    if(identify_key(s,&fearch->horizontal_diffusion2D) )  continue;
    if(identify_key(s,&fearch->horizontal_advection2D) )  continue;
    if(identify_key(s,&fearch->transport2D) )             continue;
    if(identify_key(s,&fearch->divergence2D) )            continue;
    if(identify_key(s,&fearch->u_flag) )                continue;
    if(identify_key(s,&fearch->w_flag) )                continue;
    if(identify_key(s,&fearch->T_flag) )                continue;
    if(identify_key(s,&fearch->S_flag) )                continue;
    if(identify_key(s,&fearch->R_flag) )                continue;
    if(identify_key(s,&fearch->horizontal_diffusion) )  continue;
    if(identify_key(s,&fearch->horizontal_advection) )  continue;
    if(identify_key(s,&fearch->ibd_flag) )              continue;
    if(identify_key(s,&fearch->Pa_flag) )               continue;
    if(identify_key(s,&fearch->ws_flag) )               continue;
    if(identify_key(s,&fearch->lateral_mixing) )        continue;
    if(identify_key(s,&fearch->radiation_stress) )        continue;
    if(identify_key(s,&fearch->gradSxx) )        		continue;
    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
finished:
  fearch->active=(fearch->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_structured_keyword( structured_C *structured,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    //cout << s <<flush;
    structured->active=true;
    pos=s.find("##"); if(pos!=-1)goto finished;

    if(identify_key(s,&structured->OnOffFlag) )               continue;
    if(identify_key(s,&structured->zone) )                    continue;
    if(identify_key(s,&structured->OverwriteAppendFlag) )     continue;
    if(identify_key(s,&structured->ArchiveStart) )            continue;
    if(identify_key(s,&structured->ArchiveInterval) )         continue;
    if(identify_key(s,&structured->Start) )                   continue;
    if(identify_key(s,&structured->Interval) )                continue;
    if(identify_key(s,&structured->z_flag) )                  continue;
    if(identify_key(s,&structured->ubar_flag) )               continue;
    if(identify_key(s,&structured->horizontal_diffusion2D) )  continue;
    if(identify_key(s,&structured->horizontal_advection2D) )  continue;
    if(identify_key(s,&structured->transport2D) )             continue;
    if(identify_key(s,&structured->divergence2D) )            continue;
    if(identify_key(s,&structured->zrange) )                  continue;
    if(identify_key(s,&structured->MKE2D) )                   continue;
    if(identify_key(s,&structured->Cd_flag) )                 continue;
    if(identify_key(s,&structured->z0_flag) )                 continue;
    if(identify_key(s,&structured->ibd_flag) )                continue;
    if(identify_key(s,&structured->Pa_flag) )                 continue;
    if(identify_key(s,&structured->ws_flag) )                 continue;
    if(identify_key(s,&structured->u_flag) )                  continue;
    if(identify_key(s,&structured->w_flag) )                  continue;
    if(identify_key(s,&structured->T_flag) )                  continue;
    if(identify_key(s,&structured->S_flag) )                  continue;
    if(identify_key(s,&structured->R_flag) )                  continue;
    if(identify_key(s,&structured->horizontal_diffusion) )    continue;
    if(identify_key(s,&structured->horizontal_advection) )    continue;
    if(identify_key(s,&structured->pressure_gradient) )       continue;
    if(identify_key(s,&structured->lateral_mixing) )          continue;
    if(identify_key(s,&structured->tracer2D) )                continue;
    if(identify_key(s,&structured->radiation_stress) )        continue;
   pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }
finished:
  structured->active=(structured->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_sample_keyword( sample_C *sample,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) ) {
    sample->active=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1)goto finished;

    if(identify_key(s,&sample->OnOffFlag) )      continue;
    if(identify_key(s,&sample->file) )           continue;
    if(identify_key(s,&sample->start) )          continue;
    if(identify_key(s,&sample->SaveInterval) )   continue;
    if(identify_key(s,&sample->OverAppendFlag) ) continue;

    if(identify_key(s,&sample->SectionOnOffFlag) )      continue;
    if(identify_key(s,&sample->SectionInputFile) )      continue;
    if(identify_key(s,&sample->SectionStart) )          continue;
    if(identify_key(s,&sample->SectionSaveInterval) )   continue;
    if(identify_key(s,&sample->SectionOverAppendFlag) ) continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }

finished:
  sample->active=(sample->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  sample->SectionActive=(sample->SectionOnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_analysis_keyword( analysis_C *analysis,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) ) {
    analysis->active=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1)goto finished;

    if(identify_key(s,&analysis->OnOffFlag) )         continue;
    if(identify_key(s,&analysis->ListFile) )          continue;
    if(identify_key(s,&analysis->AutoCompletion) )    continue;
    if(identify_key(s,&analysis->HarmoStart) )        continue;
    if(identify_key(s,&analysis->SampleInterval) )    continue;
    if(identify_key(s,&analysis->ComputeInterval) )   continue;
    if(identify_key(s,&analysis->Strategy) )          continue;
    if(identify_key(s,&analysis->SequentialRecycle) )          continue;
    if(identify_key(s,&analysis->SequentialArchiveLevel) )          continue;
    if(identify_key(s,&analysis->LGP1ArchiveFlag) )   continue;
    if(identify_key(s,&analysis->NETCDFArchiveFlag) ) continue;

    if(identify_key(s,&analysis->SGArchiveFlag) )     continue;
    if(identify_key(s,&analysis->grid) )              continue;
    if(identify_key(s,&analysis->z_flag) )            continue;
    if(identify_key(s,&analysis->u_flag) )            continue;
    if(identify_key(s,&analysis->LSA_flag) )          continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    }

finished:
  analysis->active=(analysis->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  analysis->SGsave=(analysis->SGArchiveFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}
//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_energy_keyword( energy_C *energy,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) ) {
    energy->active=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1)goto finished;

    if(identify_key(s,&energy->OnOffFlag) ) continue;
    if(identify_key(s,&energy->start) ) continue;
    if(identify_key(s,&energy->interval) ) continue;
    if(identify_key(s,&energy->TimeStepComp) ) continue;
    if(identify_key(s,&energy->HarmoComp) ) continue;
    if(identify_key(s,&energy->HarmoCompByRegions) ) continue;
    if(identify_key(s,&energy->ResonanceExploration) ) continue;
    if(identify_key(s,&energy->ResonanceFrequencies) ) continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;

    } 
finished:
  energy->active=(energy->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_drifters2D_keyword( drifters2D_C *drifters2D,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    drifters2D->active=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1) goto finished;

    if(identify_key(s,&drifters2D->start) ) continue;
    if(identify_key(s,&drifters2D->interval) ) continue;
    if(identify_key(s,&drifters2D->InputFile) ) continue;
    if(identify_key(s,&drifters2D->OnOffFlag) ) continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    } 
finished:
  drifters2D->active=(drifters2D->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_tracers_keyword( tracers2D_C *tracers,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    tracers->active=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1) goto finished;

    if(identify_key(s,&tracers->start) )                           continue;
    if(identify_key(s,&tracers->interval) )                        continue;
    if(identify_key(s,&tracers->TimeIntegration) )                 continue;
    if(identify_key(s,&tracers->InputFile) )                       continue;
    if(identify_key(s,&tracers->InitialisationFile) )              continue;
    if(identify_key(s,&tracers->HorizontalDiffusionMode) )         continue;
    if(identify_key(s,&tracers->HorizontalDiffusionCoefficient) )  continue;
    if(identify_key(s,&tracers->SmagorinskyCoefficient) )          continue;
    if(identify_key(s,&tracers->UpwindCoefficient) )               continue;
    if(identify_key(s,&tracers->OnOffFlag) ) continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    } 
finished:
  tracers->active=(tracers->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_tracers3D_keyword( tracers3D_C *tracers3D,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
//    tracers3D->active=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1) goto finished;

    if(identify_key(s,&tracers3D->OnOffFlag) )                 continue;
    if(identify_key(s,&tracers3D->TracersIntegration) )        continue;
    if(identify_key(s,&tracers3D->TracersMode) )               continue;
    if(identify_key(s,&tracers3D->TracersDiscretisation) )     continue;
    if(identify_key(s,&tracers3D->TracersAdvectionCentering) ) continue;
    if(identify_key(s,&tracers3D->TracerInitialisation) )      continue;
    if(identify_key(s,&tracers3D->InitialisationDirectory) )   continue;
    if(identify_key(s,&tracers3D->TFile) )                     continue;
    if(identify_key(s,&tracers3D->TVariable) )                 continue;
    if(identify_key(s,&tracers3D->SFile) )                     continue;
    if(identify_key(s,&tracers3D->SVariable) )                 continue;
    if(identify_key(s,&tracers3D->RFile) )                     continue;
    if(identify_key(s,&tracers3D->RVariable) )                 continue;
    if(identify_key(s,&tracers3D->PFile) )                     continue;
    if(identify_key(s,&tracers3D->LateralDiffusion) )          continue;
    if(identify_key(s,&tracers3D->HorizontalDiffusion) )       continue;
    if(identify_key(s,&tracers3D->VerticalDiffusion) )         continue;
    if(identify_key(s,&tracers3D->HorizontalSmooth) )          continue;
    if(identify_key(s,&tracers3D->ComputeDiagnostics) )        continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    } 
finished:
//  tracers3D->active=(tracers3D->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_streamf_keyword( streamf_C *streamf,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    streamf->active=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1)goto finished;

    if(identify_key(s,&streamf->ComputeInterval) ) continue;
    if(identify_key(s,&streamf->OnOffFlag) ) continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    } 
finished:
  streamf->active=(streamf->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 bool read_input_constraint_keyword( constraint_C *constraint,ifstream *in )

//XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int pos;
  string s;

  while(getline(*in,s) )
    {
    constraint->active=true;
    //cout << s <<flush;
    pos=s.find("##"); if(pos!=-1)goto finished;

    if(identify_key(s,&constraint->ElevationConstraintFile) ) continue;
    if(identify_key(s,&constraint->OnOffFlag) )               continue;

    pos=s.find("%"); if( (pos!=-1) ) continue;

    cout <<__FUNCTION__<< " : can not understand this line\n" << s  << "\n" << flush;continue;;
    } 
finished:
  constraint->active=(constraint->OnOffFlag.actual_value == Keywords::KEY_TRUE);
  return(true);
}



