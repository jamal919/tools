
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief tidal atlas definitions
*/
/*----------------------------------------------------------------------------*/

#include "config.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "list.h"
#include "map.h"
#include "mgr.h"
#include "grd.h"
#include "filter.h"
#include "functions.h"
#include "tides.h"
#include "poc-netcdf-data.hpp"
#include "fe-proto.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int decode_name(date_t actual, const char *varname, const char *name_template, char **filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int  status;
  char dummy[256], *pointer;
  date_t cdf_reference;
  FILE *out;

/*------------------------------------------------------------------------------
  build the meteo file name and open it*/
  (*filename) = new char[strlen(name_template) + 256];
  sprintf((*filename), "%s", name_template);

/*
  status=get_file_size(*filename, 0);
*/

  out = fopen(*filename, "r");
  status = (out == NULL);

  switch (status) {
    case 0:
/*------------------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      break;

    default:
/*------------------------------------------------------------------------------
      use format information*/
      pointer = strstr((*filename), "YYYY");
      if(pointer != NULL) {
        sprintf(dummy, "%4d", actual.year);
        strncpy(pointer, dummy, 4);
        }
      pointer = strstr((*filename), "MM");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.month);
        strncpy(pointer, dummy, 2);
        }
      pointer = strstr((*filename), "DD");
      if(pointer != NULL) {
        sprintf(dummy, "%2.2d", actual.day);
        strncpy(pointer, dummy, 2);
        }
      pointer = strstr((*filename), "VARNAME");
      if((pointer != NULL) && (varname != NULL)) {
        sprintf(dummy, "%s", varname);
        strncpy(pointer, dummy, strlen(varname));
        }
      break;
    }

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_decode_atlasname_help(int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///prints convention help for decode_atlasname()
/**
\param mode Default: 3. See decode_atlasname()
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "CONVENTION\n"
    "  \"WAVE\" is replaced by the name of the wave (which is all upper-case).\n"
    "  \"VAR\" is replaced by the name of the variable.\n");
  if(mode&2){
    printf("  \"wave\" is replaced by the lower-cased name of the wave.\n");
    printf("  \"Wave\" is replaced by the all-but-first-letter-lower-cased name of the wave.\n");
    }
  if(mode&1){
    printf("  If the file exists without any replacements made, do not do any replacement.\n");
    }
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *decode_atlasname(const char *convention,const char *wave, const char*var, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///computes a name from a convention
/**
\date 2011-11-30 Damien Allain : created

\param *convention convention : "WAVE" will be replaced with the name of the wave and "VAR" with the name of the variable
\param *wave name of the wave. If NULL no replacement will be made.
\param *var name of the variable. If NULL no replacement will be made.
\param mode Default: 3. mode&1 : check existence. mode&2 : replace "wave" with lower-cased name of the wave and "Wave" with all-but-first-letter-lower-cased name of the wave
\returns a new-allocated atlas name

\sa print_decode_atlasname_help() tide_decode_atlasname()
*/
/*----------------------------------------------------------------------------*/
{
  FILE *out;
  if(mode&1 && (out = fopen(convention, "r"))!=NULL){
    fclose(out);
    return poc_strdup(convention);
    }
  
  string b(convention);//buffer
  const int nC=4;
  const char *ns[nC]={"wave","Wave","WAVE","VAR"};//needles
  char *r;//replacement string
  int i,nI,nL,rL;//char index, needle index, needle length, replacement length
  
  for(nI=0;nI<nC;nI++){
    switch(nI){
      case 0:
      case 1:
        if(wave==0)continue;
        if(!(mode&2))
          continue;
        r=poc_strdup(wave);
        for(i=nI;r[i];i++)r[i]=tolower(r[i]);
        break;
      case 2:
        if(wave==0)continue;
        r=poc_strdup(wave);break;
      case 3:
        if(var==0)continue;
        r=poc_strdup(var);break;
      }
    rL=strlen(r);
    if(rL==0)goto delete_r;
    nL=strlen(ns[nI]);
    while((i=b.find(ns[nI]))>=0){
      b.replace(i,nL,r,rL);
      }
    delete_r:
    delete[]r;
    }
  
  return poc_strdup(b.c_str());
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_tide_decode_atlasname_help(int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///prints convention help for tide_decode_atlasname()
/**
\param mode Default: 3. See tide_decode_atlasname()
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "CONVENTION\n"
    "  \"WAVE\" is replaced by the name of the wave (which is all upper-case).\n");
  if(mode&2){
    printf("  \"wave\" is replaced by the lower-cased name of the wave.\n");
    printf("  \"Wave\" is replaced by the all-but-first-letter-lower-cased name of the wave.\n");
    }
  if(mode&1){
    printf("  If the file exists without any replacements made, do not do any replacement.\n");
    }
  /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_decode_atlasname(const char *atlas_directory, const char *atlas_convention,
                            const char *wave, const char *varname, char **filename, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///computes a name from a convention
/**
\param *atlas_directory directory convention : see atlas_convention below for the replacements
\param *atlas_convention file name convention : "WAVE" will be replaced with the name of the wave (which is all upper-case, but see mode below)
\param *wave name of the wave
\param **filename pointer to computed file name. Will be allocated with \b new
\param mode Default: 3.
mode&1 : if the file exists without any replacements made, do not do any replacement and return 0.
mode&2 : replace "wave" with lower-cased name of the wave and "Wave" with all-but-first-letter-lower-cased name of the wave
\return 0 if file is found or -1

\sa print_tide_decode_atlasname_help() decode_atlasname()
*/
/*----------------------------------------------------------------------------*/
{
  int i;
  char dummy[256], *pointer=NULL, *p2=NULL;
  FILE *out=NULL;

  if(atlas_convention==0) {
    STDERR_BASE_LINE("atlas naming template not given, quit \n");
    return(-1);
    }
/*------------------------------------------------------------------------------
  build the atlas file name*/
  *filename=0;
  if(atlas_directory==0) {
    (*filename) = new char[strlen(atlas_convention)+10];
    sprintf((*filename), "%s", atlas_convention);
    }
  else {
    (*filename) = new char[strlen(atlas_directory) + 1 + strlen(atlas_convention)+10];
    sprintf((*filename), "%s/%s", atlas_directory, atlas_convention);
    }

  if(mode&1){
/*------------------------------------------------------------------------------
    check existence name*/
    out = fopen(*filename, "r");
    if(out!=NULL){
/*------------------------------------------------------------------------------
      file exists, do nothing more*/
      fclose(out);
      return 0;
      }
    }

  ///\note Only the first occurrence of WAVE, wave or Wave, in this order, is replaced!
  for(i=0;i<3;i++){
/*------------------------------------------------------------------------------
    use format information*/
    switch(i){
    case 0:
      pointer = strstr((*filename), "WAVE");
      break;
    case 1:
      pointer = strstr((*filename), "wave");
      break;
    case 2:
      pointer = strstr((*filename), "Wave");
      break;
      }
    if(pointer != NULL) {
      p2 =strdup(pointer);
      sprintf(dummy, "%s", wave);
      int k,kMax=strlen(dummy);
      switch(i){
      case 0:
        k=kMax;
        break;
      case 1:
        k=0;
        break;
      case 2:
        k=1;
        break;
        }
      for(;k<kMax;k++) {
        dummy[k]=tolower(dummy[k]);
        }
      strcpy(pointer, dummy);
      strcat(pointer, p2+4);
      free(p2);
      break;
      }
    if(!(mode&2))
      break;
    }
  
  if(varname!=0) {
    pointer = strstr((*filename), "VARNAME");
    if(pointer != NULL) {
      p2 =strdup(pointer);
      sprintf(dummy, "%s", varname);
      strncpy(pointer, dummy, strlen(varname));
      strcat(pointer, p2+7);
      }
    }

  /* check existence */
  if ((out = fopen(*filename, "r"))) {
    fclose(out);
    return 0;
  }
  return -1;
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int tide_decode_varname(const char *vtemplate, const char *wave, char **varname)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int i;
//   char dummy[256], *pointer=NULL, *p2=NULL;
//   FILE *out=NULL;
// 
//   if(vtemplate==0) {
//     STDERR_BASE_LINE("variable naming template not given, quit \n");
//     return(-1);
//     }
//     
//   (*varname) = new char[strlen(vtemplate)+10];
//   sprintf((*varname), "%s", vtemplate);
// 
//   ///\note Only the first occurrence of WAVE, wave or Wave, in this order, is replaced!
//   for(i=0;i<3;i++){
// /*------------------------------------------------------------------------------
//     use format information*/
//     switch(i){
//       case 0:
//         pointer = strstr((*varname), "WAVE");
//         break;
//       case 1:
//         pointer = strstr((*varname), "wave");
//         break;
//       case 2:
//         pointer = strstr((*varname), "Wave");
//         break;
//         }
//       
//     if(pointer != NULL) {
//       p2 =strdup(pointer);
//       sprintf(dummy, "%s", wave);
//       int k,kMax=strlen(dummy);
//       switch(i){
//       case 0:
//         k=kMax;
//         break;
//       case 1:
//         k=0;
//         break;
//       case 2:
//         k=1;
//         break;
//         }
//       for(;k<kMax;k++) {
//         dummy[k]=tolower(dummy[k]);
//         }
//       strcpy(pointer, dummy);
//       strcat(pointer, p2+4);
//       free(p2);
//       break;
//       }
// //     if(!(mode&2))
// //       break;
//     }
//    
//   return(0);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int tide_atlas2mgr_template(const char *filename, const char *v1, const char *v2, vector<mgr_t> mgr, int nmgr, T *a,T *G, T mask,atlas_grid_or_mesh_t *gm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
///extract amplitudes and phases of structured atlas for the given tide gauges
{
  int n,status;
  double *lat,*lon;//lat and lon of the given locations
  lat=new double[nmgr];
  lon=new double[nmgr];

  for(n = 0; n < nmgr; n++) {
    if(strcmp(mgr[n].loc.units,"degrees")==0) {
      lon[n] = mgr[n].loc.lon;
      lat[n] = mgr[n].loc.lat;
      }
    else {
      lon[n] = mgr[n].loc.lon * r2d;
      lat[n] = mgr[n].loc.lat * r2d;
      }
    }

  status=tide_atlas2positions(filename, v1, v2, lon, lat, nmgr, a, G, mask,gm);
  delete[]lat;
  delete[]lon;
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_atlas2mgr(const char *filename, const char *v1, const char *v2, vector<mgr_t> mgr, int nmgr, float *a,float *G, float mask,atlas_grid_or_mesh_t *gm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  extract amplitudes and phases of structured atlas at given locations
-----------------------------------------------------------------------------*/
{
  int status;
  status=tide_atlas2mgr_template(filename, v1, v2, mgr, nmgr, a, G, mask, gm);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_atlas2mgr(const char *filename, const char *v1, const char *v2, vector<mgr_t> mgr, int nmgr, double *a,double *G, double mask,atlas_grid_or_mesh_t *gm)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  extract amplitudes and phases of structured atlas at given locations
-----------------------------------------------------------------------------*/
{
  int status;
  status=tide_atlas2mgr_template(filename, v1, v2, mgr, nmgr, a, G, mask, gm);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int tide_atlasSG2positions_template(const char *filename, const char *v1, const char *v2,const double *lon,const double *lat, int npositions, T *a,T *G, T mask,const grid_t &grid,int verbose,const char *wave, int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  extract amplitudes and phases of structured atlas at given locations
-----------------------------------------------------------------------------*/
{
  complex<T> cmask,*z=NULL;
  int status,n;

  cmask=complex<T>(9999.,9999.);
  z=new complex<T>[npositions];

  status=tide_atlasSG2positions(filename, v1, v2, lon, lat, npositions, z, cmask,grid,verbose,wave,strict);
  if(status!=NC_NOERR) goto finished;

  for(n = 0; n < npositions; n++) {
    if(z[n]==cmask) {
      a[n] =  mask;
      G[n] =  mask;
      }
    else {
      a[n] =  abs(z[n]);
/*------------------------------------------------------------------------------
      convert into phase lag : 0 < G < 360 */
      G[n] = -arg(z[n]) * r2d;
      while(G[n] < 0){
        G[n] += 360;
        }
      while(G[n] > 360.0){
        G[n] -= 360;
        }
      }
    }

finished:
  delete[]z;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_atlasSG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, float *a,float *G, float mask,const grid_t &grid,int verbose,const char *wave,int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlasSG2positions_template(filename, v1, v2, lon, lat, npositions, a,G, mask,grid,verbose,wave,strict);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_atlasSG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, double *a,double *G, double mask,const grid_t &grid,int verbose,const char *wave,int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlasSG2positions_template(filename, v1, v2, lon, lat, npositions, a,G, mask,grid,verbose,wave,strict);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int tide_atlas2positions_template(const char *filename, const char *v1, const char *v2,const double *lon,const double *lat, int npositions, T *a,T *G, T mask,atlas_grid_or_mesh_t *gm,int verbose,const char *wave, int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  extract amplitudes and phases of structured atlas at given locations
-----------------------------------------------------------------------------*/
{
  complex<T> cmask,*z=NULL;
  int status,n;

  cmask=complex<T>(9999.,9999.);
  z=new complex<T>[npositions];

  status=tide_atlas2positions(filename, v1, v2, lon, lat, npositions, z, cmask,gm,verbose,wave,strict);

  for(n = 0; n < npositions; n++) {
    if(z[n]==cmask) {
      a[n] =  mask;
      G[n] =  mask;
      }
    else {
      a[n] =  abs(z[n]);
/*------------------------------------------------------------------------------
      convert into phase lag : 0 < G < 360 */
      G[n] = -arg(z[n]) * r2d;
      while(G[n] < 0){
        G[n] += 360;
        }
      while(G[n] > 360.0){
        G[n] -= 360;
        }
      }
    }

  delete[]z;
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_atlas2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, float *a,float *G, float mask,atlas_grid_or_mesh_t *gm,int verbose,const char *wave,int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlas2positions_template(filename, v1, v2, lon, lat, npositions, a,G, mask,gm,verbose,wave,strict);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tide_atlas2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, double *a,double *G, double mask,atlas_grid_or_mesh_t *gm,int verbose,const char *wave,int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlas2positions_template(filename, v1, v2, lon, lat, npositions, a,G, mask,gm,verbose,wave,strict);
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

hconstant_t *tide_atlas2positions(const char *atlPathConv,const spectrum_t &WaveList,const char *ampNameConv,const char *phaNameConv,double *lon,double *lat,int n,const double mask,int verbose,int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// looped wrapper for tide_atlas2position(const char*,const char*,const char *,double*,double*,int,double*a,double*G,double,int)
/**
\param *atlPathConv atlas path convention
\param *WaveList
\param *ampNameConv amplitude variable name convention
\param *phaNameConv phase variable name convention
\param *lon
\param *lat
\param *n number of points
\param mask
\param verbose \c verbose&3 is given as verbosity to tide_atlas2positions(), if \c verbose&4, prints constants
\returns a hconstant_t[n] array, initialised with hconstant_t::init_polar()
*/
/*----------------------------------------------------------------------------*/
{
  hconstant_t *constants;
  double *a,*G;//<amplitudes and phases in the atlas at the points
  char *atlPath,*ampName,*phaName;//<decoded names of atlas,amplitudes and phases
  atlas_grid_or_mesh_t gm;
  char *wavename;
  int i,j;//wave and point indexes
  
  int constant_verbosity,tide_atlas2positions_verbosity;
  constant_verbosity=verbose&4;
  tide_atlas2positions_verbosity=verbose&3;
  verbose=constant_verbosity and tide_atlas2positions_verbosity>=0;
  
  constants=new hconstant_t[n];
  for(j=0;j<n;j++){
    constants[j].init_polar(WaveList.n);
    }
  
  a=new double[n];
  G=new double[n];
  
  for(i=0;i<WaveList.n;i++){
    wavename=WaveList.waves[i].name;
    
    tide_decode_atlasname(NULL,atlPathConv,wavename, 0, &atlPath,3);
    tide_decode_atlasname(NULL,ampNameConv,wavename, 0, &ampName,2);
    tide_decode_atlasname(NULL,phaNameConv,wavename, 0, &phaName,2);
    
    if(verbose)printf("(%d/%d)Reading constants for %10s (%8.5g deg/h) from %s :\n",i+1,WaveList.n,wavename,WaveList.waves[i].omega,atlPath);fflush(stdout);
    
    //read constants in the atlases
    tide_atlas2positions(atlPath, ampName, phaName, lon, lat, n, a,G, NAN, &gm, tide_atlas2positions_verbosity, wavename,strict);
    
    if(constant_verbosity)
      printf("Wave : %s\npoint|longitude|latitude|amplitude|phase(deg)|complex\n",wavename);
    
    for(j=0;j<n;j++){
      complex<double> z=polar(a[j],G[j]*d2r);
      constants[j].a[i]=a[j];
      constants[j].G[i]=G[j];
      
      if(constant_verbosity)printf("%5d|%9.3g|%8.3g|%9.3g|%-+10.4g|%g%+gj\n",j,lon[j],lat[j],constants[j].a[i],constants[j].G[i],real(z),imag(z));
      }
    
    if(constant_verbosity)fflush(stdout);
    
    delete[]atlPath;
    delete[]ampName;
    delete[]phaName;
    }
  
  delete[]a;
  delete[]G;
  
  return constants;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void scale_constants(hconstant_t *constants,int pn,int wn,double factor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// mutiply constants by a factor
/**
\param[in,out] *constants a hconstant_t[n] array, initialised with hconstant_t::init_polar()
\param pn number of points
\param wn number of waves
\param factor
*/
/*----------------------------------------------------------------------------*/
{
  int pI,wI;
  
  if(factor==1.)
    return;
  
  for(pI=0;pI<pn;pI++){
    hconstant_t *constanti=&constants[pI];
    for(wI=0;wI<wn;wI++)
      constanti->a[wI]*=factor;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void tide_atlas2positions_fortran_template(const char *aPC,int *aPCl,spectrum_t *WaveList,const char *avC,int *avl,const char *gvC,int *gvl,const double *lon,const double *lat, int *npositions, T *a,T *G, T *mask,int *verbose,int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran interpolator template
/**
\param *aPC atlas path convention
\param *aPCl length of above
\param *WaveList
\param *avC amplitude variable name convention
\param *avl length of above
\param *gvC phase variable name convention
\param *gvl length of above
\param *lon
\param *lat
\param *npositions
\param *a amplitudes
\param *G phases
\param *mask
\param *verbose
\param *status

\note \c aPC, \c avC and \c agC MUST BE SPACE-TERMINATED STRINGS BECAUSE OF C TO FORTRAN MESS !
\note actually Fortran strings are space padded anyway :)
*/
/*----------------------------------------------------------------------------*/
{
  int status_;
  if(status==0){
    STDOUT_BASE_LINE("WARNING:STATUS POINTER IS %p. PATCHING...\n",status);
    status=&status_;
    }
  
  /* sort out SPACE TERMINATED Fortran strings */
  const string
    aPC_(aPC,*aPCl),
    avC_(avC,*avl),
    gvC_(gvC,*gvl);
  
  char *aP,*av,*gv;
  char *wavename;
  int i,m;//wave index, constant index
  
  atlas_grid_or_mesh_t gm;
  
  if(*verbose){
    range_t<double> loR;
    const double lonErr=36E3;
    
    STDOUT_BASE_LINE("testing lon...");
    for(i=0;i<*npositions;i++){
      const double loni=lon[i];
      if(loni<-lonErr or loni>lonErr) TRAP_ERR_EXIT(-1,"lon[%d]=%g\n",i,loni);
      loR<<loni;
      }
    printf(" in [%g;%g]\n",loR.min,loR.max);
    }
  
  for(i=0;i<WaveList->n;i++){
    wavename=WaveList->waves[i].name;
    
    m=i* *npositions;
    
    ///decode conventions with tide_decode_atlasname()
    tide_decode_atlasname(NULL,aPC_.c_str(),wavename, 0, &aP);
    tide_decode_atlasname(NULL,avC_.c_str(),wavename, 0, &av);
    tide_decode_atlasname(NULL,gvC_.c_str(),wavename, 0, &gv);
    
    if(*verbose) STDOUT_BASE_LINE("Reading constants for %10s (%8.5g deg/h) from %s(%s,%s)\n",wavename,WaveList->waves[i].omega,aP,av,gv);
    ///read the constants in the atlases with tide_atlas2positions_template()
    *status=tide_atlas2positions_template(aP, av,gv, lon, lat, *npositions, &a[m],&G[m], *mask, &gm, *verbose, wavename, 0);
    
    if(*status!=0) NC_CHKERR_BASE_LINE(*status,"tide_atlas2positions_template(\"%s\",\"%s\",\"%s\",,,,,,,,\"%s\",0) error",aP,av,gv,wavename);
    
    delete[]aP;
    delete[]av;
    delete[]gv;
    
    if(*status!=0) break;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void tide_atlas2positions4_(const char *aPC,int *aPCl,spectrum_t *WaveList,const char *avC,int *avl,const char *gvC,int *gvl,const double *lon,const double *lat,int *npositions,float *a,float *G, float *mask,int *verbose,int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran interpolator
/*----------------------------------------------------------------------------*/
{
  tide_atlas2positions_fortran_template(aPC,aPCl,WaveList,avC,avl,gvC,gvl,lon,lat,npositions,a,G,mask,verbose,status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void tide_atlas2positions8_(const char *aPC,int *aPCl,spectrum_t *WaveList,const char *avC,int *avl,const char *gvC,int *gvl,const double *lon,const double *lat,int *npositions,double *a,double *G, double *mask,int *verbose,int *status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran interpolator
/*----------------------------------------------------------------------------*/
{
  tide_atlas2positions_fortran_template(aPC,aPCl,WaveList,avC,avl,gvC,gvl,lon,lat,npositions,a,G,mask,verbose,status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int tide_atlasSG2positions_template(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions,complex<T> *z,complex<T> cmask,const grid_t &grid,int verbose,const char *wave, int strict)
//similar to extract_SGatlas
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///extract complex amplitudes of structured atlas at given locations
/**
\param verbose Default:1
\param strict disable extrapolation
\sa tide_atlas2positions_template() tide_atlasUG2positions_template()
*/
/*----------------------------------------------------------------------------*/
{
  int k, l, n, status;
  int kk,ll,mm;
  const int nnext=24;
  int nextk[nnext]={-1,+1, 0, 0,-1,+1,-1,+1,-2,+2, 0, 0,-2,+2,-2,+2,-1,-1,+1,+1,-2,-2,+2,+2};
  int nextl[nnext]={ 0, 0,-1,+1,-1,+1,+1,-1, 0, 0,-2,+2,-1,-1,+1,+1,-2,+2,-2,+2,-2,+2,-2,+2};
  int64_t accel=-1;
  complex<T> zz;
  complex<T> *tide=NULL;
  double x, y;
  int frame=-1;
  int nerrors=0;

  const int gridHsize=grid.Hsize();
  exitIfNull(tide   =new complex<T>[gridHsize]);
  
  if(wave!=NULL){
    status=poc_get_frame_from_name(filename,wave,&frame,verbose);
    switch(status){
      case 0:break;
      case -1:NC_TRAP_ERROR(return,NC_EINVALCOORDS,verbose,"poc_get_frame_from_name(\"%s\",\"%s\",) returned %d, so returning",filename,wave,status);
      default:NC_TRAP_ERROR(return,NC_ENOTVAR,verbose,"poc_get_frame_from_name(\"%s\",\"%s\",) returned %d, so returning",filename,wave,status);
      }
    }

/*------------------------------------------------------------------------------
  load netcdf variable */
  status=poc_get_cvara(filename,v1,v2,frame,tide,verbose);
  const complex<T> cmask0=NC_FILL_DOUBLE;
  for(n=0;n<gridHsize;n++){
    complex<T> *tiden=&tide[n];
    if(*tiden==cmask0)
      *tiden=cmask;
    }
  
/*------------------------------------------------------------------------------
  interpolate */
  
  for(n = 0; n < npositions; n++) {
    x = lon[n];
    y = lat[n];
//     if(x < grid.xmin - grid.dx / 2.)
//       x = x + 360.0;
//     if(x > grid.xmax + grid.dx / 2.)
//       x = x - 360.0;
    if(x < grid.xmin)
      x = x + 360.0;
    if(x > grid.xmax)
      x = x - 360.0;
    if(strict!=0){
      index_interpolation(grid,x,y,&accel, tide, cmask,&zz);
      }
    else{
      status = map_interpolation(grid, tide, cmask, x, y,&zz);
      if((status != 0) || (zz == cmask)) {
        status=map_index( grid,  x,  y, &k, &l);
        if(status == 0) {
          int kn;
          status=-1;
          for(kn=1;kn<=2;kn++){
          for(mm=0;((status != 0) || (zz == cmask))&&(mm<nnext);mm++){
            kk=k+kn*nextk[mm];
            ll=l+kn*nextl[mm];
            if(kk>=grid.nx || ll>=grid.ny || kk<0 || ll<0)continue;
            x=map_grid_x(grid,kk,ll);
            y=map_grid_y(grid,kk,ll);
            status = map_interpolation(grid, tide, cmask, x, y,&zz);
            }
            }
          }
        }
      }
    if(zz == cmask) {
//      printf("interpolation error: node %d lon=%lf lat=%lf \n", n, x, y);
      z[n] =  cmask;
      nerrors++;
      }
    else {
      z[n] =  zz;
      }
    }

  delete[]tide;

  status=0;
  STDOUT_BASE_LINE("%d points = %d masked + %d interpolated \n", npositions, nerrors, npositions-nerrors);
  fflush(stdout);
  return (status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int tide_atlasSG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions,complex<double> *z,complex<double> cmask,const grid_t &grid,int verbose,const char *wave,int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlasSG2positions_template(filename,v1,v2,lon,lat,npositions,z,cmask,grid, verbose, wave, strict);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int tide_atlasSG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions,complex<float> *z,complex<float> cmask,const grid_t &grid,int verbose,const char *wave, int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlasSG2positions_template(filename,v1,v2,lon,lat,npositions,z,cmask,grid, verbose, wave, strict);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int tide_atlas2positions_template(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions,complex<T> *z,complex<T> cmask,atlas_grid_or_mesh_t *gm,int verbose,const char *wave, int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///extract complex amplitudes of structured or unstructured atlas at given locations
/**
\param verbose Default:1
\param strict for tide_atlasSG2positions_template()
calls tide_atlasSG2positions_template() then tide_atlasUG2positions_template()
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  poc_grid_data_t gdata;
  
  poc_global_t global;
  poc_var_t var;
  
  status=poc_inq_var(filename,v1,&var);
  
/*------------------------------------------------------------------------------
  check grid or mesh */
  struct timeval before;
  gettimeofday(&before);
  
  if(verbose>0)STDERR_BASE_LINE("poc_get_grid_data(\"%s\",(\"%s\"),...)\n",filename,v1);
  
  status=poc_get_grid_data(filename,var,&gdata,verbose-1,&global);
  if(status)NC_TRAP_ERROR(return,status,verbose,"poc_get_grid_data(\"%s\",...) error",filename);
  
  if(gdata!=gm->data){
    /* update grid or mesh */
    bool isGrid,isMesh;
    
    gm->destroy();
    
    gm->data.transfer(&gdata);
    
    status=test_grid_data(gm->data,&isGrid,&isMesh,verbose-1);
    NC_TRAP_ERROR(return,status,verbose,"test_grid_data() error");
    
    if(isGrid){
      if(verbose>0)
        STDERR_BASE_LINE("%gs:poc_grid_data_to_grid(,,,[\""+gm->data.xv.info.name+"\",\""+gdata.yv.info.name+"\"],...)\n",difftime(before));
      
      status=poc_grid_data_to_grid(filename,global,var,gm->data,&gm->grid,verbose-1,-1);
      
      if(strict!=0){
        set_grid_list(&gm->grid, verbose);
        }
      
      }
    else if(isMesh){
      if(verbose>0)
        STDERR_BASE_LINE("%gs:poc_grid_data_to_mesh ...\n",difftime(before));
      
      status=poc_grid_data_to_mesh(filename,global,var,gm->data,&gm->mesh,verbose-1,-1);
      }
    else
      TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
    }
  
/*------------------------------------------------------------------------------
  interpolate */
  
  if(gm->mesh.vertices!=0){
    gm->discretisation=poc_get_discretisation(var,&gm->mesh,verbose);
    discretisation_init(&gm->mesh,gm->discretisation,0);
    
    status=tide_atlasUG2positions(filename,v1,v2,lon,lat,npositions,z,cmask,gm->mesh,gm->discretisation,verbose-1,wave);
    }
  else{
    status=tide_atlasSG2positions_template(filename,v1,v2,lon,lat,npositions,z,cmask,gm->grid,verbose-1,wave,strict);
    }
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int tide_atlas2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat,int npositions,complex<double> *z,complex<double> cmask,atlas_grid_or_mesh_t *gm,int verbose,const char *wave,int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlas2positions_template(filename,v1,v2,lon,lat,npositions,z,cmask,gm,verbose,wave,strict);
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int tide_atlas2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat,int npositions,complex<float> *z,complex<float> cmask,atlas_grid_or_mesh_t *gm,int verbose,const char *wave,int strict)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  status=tide_atlas2positions_template(filename,v1,v2,lon,lat,npositions,z,cmask,gm,verbose,wave,strict);
  return status;
}


hconstant_t * atlas_constantsN(const char *atlas_convention,const spectrum_t &s, double *x, double *y, int nb_pos);


/*----------------------------------------------------------------------------*/
/// tidal atlas constants at a given position using admittance
/**
\date last reviewed 3 Aug 2011
\author Damien Allain

\param *atlas_directory path to atlas directory, see tide_decode_atlasname()
\param *atlas_convention atlas convention, see tide_decode_atlasname()
\param s spectrum
\param x longitude
\param y latitude
\returns the elevation constants
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  hconstant_t atlas_constants(const char *atlas_convention,const spectrum_t &s, double x, double y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  hconstant_t *elevation;
  
  elevation=atlas_constantsN(atlas_convention,s,&x,&y,1);

  return(elevation[0]);
}


/*----------------------------------------------------------------------------*/
/// tidal atlas constants at a given position
/**
\author Damien Allain
\author Sara Fleury

\param *atlas_convention atlas convention, see tide_decode_atlasname()
\param s spectrum
\param x *longitude
\param y *latitude
\param nb_pos      number of positions to consider
\returns the elevation constants
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

hconstant_t * atlas_constantsN(const char *atlas_convention,const spectrum_t &s, double *x, double *y, int nb_pos)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,status;
  int *available = aset(s.n,1);
  char *atlas_file;
  const char *ampVar=0,*phaVar=0;
  float mask=9999;  /* !!! HERE !!! was not init. Is 9999 OK ? - sara */
  hconstant_t *elevation;
  int n;

  for (k=0;k<s.n;k++) {
/* *-----------------------------------------------------------------------------
    read tidal atlas database */
    if (!strlen( s.waves[k].name)) {
      printf ("No wave name for wave No %d ??\n", k);
      available[k]=0;
      continue;
      }
    status=tide_decode_atlasname(NULL,atlas_convention,s.waves[k].name, 0,&atlas_file);
    if (status==-1) {
      available[k]=0;
      continue;
      }
    
    if(ampVar==0 or phaVar==0){
      poc_global_t glob;
      const int nVarList=3;
      const char *ampVarList[nVarList]={"Ha","elevation_a","amplitude"},*phaVarList[nVarList]={"Hg","elevation_G","phase"};
      int i;
      
      status=poc_inq(atlas_file,&glob);
      
      for(i=0;i<nVarList;i++){
        status=glob.variables.find(ampVarList[i]);
        if(status>=0) break;
        }
      if(i>=nVarList) TRAP_ERR_EXIT(ENOEXEC,"Please complete ampVarList\n");
      ampVar=ampVarList[i];
      
      for(i=0;i<nVarList;i++){
        status=glob.variables.find(phaVarList[i]);
        if(status>=0) break;
        }
      if(i>=nVarList) TRAP_ERR_EXIT(ENOEXEC,"Please complete phaVarList\n");
      phaVar=phaVarList[i];
      }
    
    STDOUT_BASE_LINE("tides01: treating %s wave, from %s (%s,%s)\n",s.waves[k].name,atlas_file,ampVar,phaVar);
    }
  
  elevation=tide_atlas2positions(atlas_convention,s,ampVar,phaVar,x,y,nb_pos,mask,0,0);

  for (n=0; n<nb_pos; n++) {

/**----------------------------------------------------------------------------
    infer missing constituents using admmittance*/
    admittance_t admittance;

    delete[] admittance_mts_check(&admittance, s, available);
    admittance_mts_compute(admittance, s, elevation[n], available, 0.0, false, 0);

    for (k=0;k<s.n;k++) {
/**----------------------------------------------------------------------------
      special case for Sa and Ssa: infer missing constituents using equilibrium*/
      if( (available[k] == 0) && (strcmp(s.waves[k].name,"Sa") == 0) ) {
        status=tidal_equilibrium(s.waves[k],x[n],y[n],&(elevation[n].a[k]),&(elevation[n].G[k]));
        }
      if( (available[k] == 0) && (strcmp(s.waves[k].name,"Ssa") == 0) ) {
        status=tidal_equilibrium(s.waves[k],x[n],y[n],&(elevation[n].a[k]),&(elevation[n].G[k]));
        }
      }
    
    admittance_mts_terminate(admittance);
    }

  return(elevation);
}

