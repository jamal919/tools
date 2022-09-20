
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada
\author  Sara Fleury        LEGOS/CNRS, Toulouse, France

\brief stugo fortran plug-in
*/
/*----------------------------------------------------------------------------*/

#include "config.h"

#include <stdio.h>
#include <unistd.h> //access
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "functions.h"
#include "poc-netcdf.hpp"

#include "flags.def"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void spectral_tugo_fortran_template_(const char *paramP,int *paramPl,const char *deltaP,int *deltaPl,const char *outD,int *outDl,const spectrum_t *wl,T *ha,T *hg,T *ua,T *ug,T *va,T *vg,int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran stugo template
/**
\param *paramP parameter file path
\param *deltaP delta parameter file path
\param *outD directory where the atlases will be written
\param *wl list of waves to read (more can be generated)
\param *ha elevation amplitudes
\param *hg elevation phases
\param *ua X velocity amplitudes
\param *ug X velocity phases
\param *va Y velocity amplitudes
\param *vg Y velocity phases

mask value is NC_FILL_DOUBLE+0i
\note \c paramP and \c deltaP MUST BE SPACE-TERMINATED STRINGS BECAUSE OF C TO FORTRAN MESS !
\note actually Fortran strings are space padded anyway :)
*/
/*----------------------------------------------------------------------------*/
{
  /* sort out SPACE TERMINATED Fortran strings */
  const string
    paramP_(paramP,*paramPl),
    deltaP_(deltaP,*deltaPl),
    outD_(outD,*outDl);
  size_t hlength=-1,ulength=-1,vlength=-1;
  
  if(paramP_!=""){
    /*--------------------------------------------------------------------------
      compose command */
      string cmd=TOOLS_BUILDDIR "/stugo";
      
      if(access(cmd.c_str(),X_OK)!=0){
        cmd=TOOLS_BINDIR "/stugo";
        }
      
      cmd+=" "+paramP_;
      
      if(deltaP_!=""){
        cmd+=" -delta "+deltaP_;
        }
      
      cmd+=" --output-path "+outD_;
      
    /*--------------------------------------------------------------------------
      run command */
      STDERR_BASE_LINE("system(\""+cmd+"\") ...\n");
#if GCC_VERSION >= 40500  /* Test for GCC >=4.5.0 */
#pragma GCC diagnostic push
/* __attribute__((unused)) seams to fail for function declarations */
#pragma GCC diagnostic ignored "-Wunused-value"
#endif
      *status=system(cmd.c_str());
      if(*status)NC_TRAP_ERROR(return;,*status,1,"system(\""+cmd+"\") error");
      }
  
/*------------------------------------------------------------------------------
  read output */
  for(int i=0;i<wl->n;i++){
    const string
      waveName=wl->waves[i].name,
      fileP=outD_+"/"+waveName+".spectral-CGRID.nc";
    string aName,gName;
    
    STDERR_BASE_LINE("reading "+fileP+"\n");
    
    aName=waveName+"_elevation_a";
    gName=waveName+"_elevation_G";
    if(hlength<0){
      *status=poc_get_var_length(fileP,aName,&hlength);
      if(*status)NC_TRAP_ERROR(return;,*status,1,"poc_get_var_length(\""+fileP+"\",\""+aName+"\",) error");
      }
    *status=poc_get_var(fileP,aName,&ha[i*hlength]);
    if(*status)NC_TRAP_ERROR(return;,*status,1,"poc_get_var(\""+fileP+"\",\""+aName+"\",) error");
    *status=poc_get_var(fileP,gName,&hg[i*hlength]);
    if(*status)NC_TRAP_ERROR(return;,*status,1,"poc_get_var(\""+fileP+"\",\""+gName+"\",) error");
    
    aName=waveName+"_u_a_u";
    gName=waveName+"_u_G_u";
    if(ulength<0){
      *status=poc_get_var_length(fileP,aName,&ulength);
      if(*status)NC_TRAP_ERROR(return;,*status,1,"poc_get_var_length(\""+fileP+"\",\""+aName+"\",) error");
      }
    *status=poc_get_var(fileP,aName,&ua[i*ulength]);
    if(*status)NC_TRAP_ERROR(return;,*status,1,"poc_get_var(\""+fileP+"\",\""+aName+"\",) error");
    *status=poc_get_var(fileP,gName,&ug[i*ulength]);
    if(*status)NC_TRAP_ERROR(return;,*status,1,"poc_get_var(\""+fileP+"\",\""+gName+"\",) error");
    
    aName=waveName+"_v_a_v";
    gName=waveName+"_v_G_v";
    if(vlength<0){
      *status=poc_get_var_length(fileP,aName,&vlength);
      if(*status)NC_TRAP_ERROR(return;,*status,1,"poc_get_var_length(\""+fileP+"\",\""+aName+"\",) error");
      }
    *status=poc_get_var(fileP,aName,&va[i*vlength]);
    if(*status)NC_TRAP_ERROR(return;,*status,1,"poc_get_var(\""+fileP+"\",\""+aName+"\",) error");
    *status=poc_get_var(fileP,gName,&vg[i*vlength]);
    if(*status)NC_TRAP_ERROR(return;,*status,1,"poc_get_var(\""+fileP+"\",\""+gName+"\",) error");
#if GCC_VERSION >= 40500  /* Test for GCC >=4.5.0 */
#pragma GCC diagnostic pop
#endif
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void spectral_tugo4_(const char *paramP,int *paramPl,const char *deltaP,int *deltaPl,const char *outD,int *outDl,const spectrum_t *wl,float *ha,float *hg,float *ua,float *ug,float *va,float *vg,int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran stugo plug-in
/*----------------------------------------------------------------------------*/
{
  ///\sa spectral_tugo_fortran_template_()
  spectral_tugo_fortran_template_(paramP,paramPl,deltaP,deltaPl,outD,outDl,wl,ha,hg,ua,ug,va,vg,status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void spectral_tugo8_(const char *paramP,int *paramPl,const char *deltaP,int *deltaPl,const char *outD,int *outDl,const spectrum_t *wl,double *ha,double *hg,double *ua,double *ug,double *va,double *vg,int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran stugo plug-in
/*----------------------------------------------------------------------------*/
{
  ///\sa spectral_tugo_fortran_template_()
  spectral_tugo_fortran_template_(paramP,paramPl,deltaP,deltaPl,outD,outDl,wl,ha,hg,ua,ug,va,vg,status);
}
