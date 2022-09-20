
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief miscallenous tidal definitions
*/
/*----------------------------------------------------------------------------*/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <iostream>

#include "tools-define.h"
#include "tools-structures.h"

#include "poc-time.h"
#include "statistic.h"

#include "tides.def"

#include "spectrum.h"

#ifndef VERBOSE
#define VERBOSE -1
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void astronomic_angle(astro_angles_t *astro_angles, double tj, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief Initialises an astro_angles_t

\param tj time elapsed between 1900-01-01 00:00 and the reference in julian centuries (of 36525 days)
  \note 1900-01-01 00:00  is the time reference for astronomical abaques
*/
/*----------------------------------------------------------------------------*/
{
  double days;
  double cosI,p,q;
  double t2,t4,sin2I,s2;
  double tgI2,P;
  double sh_tgn2,at1,at2;

  astro_angles->t_julian=tj;
  
/** From Schureman (1940, Fig1 & Tab1) #WAVELISTURL.
Longitudes are with reference to the vernal (i.e. spring) equinox,
i.e. the intersection of the Ecliptic (i.e. plane of Earth's orbit) and the Equator.
*/
  
/*------------------------------------------------------------------------------
  Sh_n Longitude of ascending lunar node */

  astro_angles->sh_N=(259.1560563 -omega_n *tj)*d2r;

/*------------------------------------------------------------------------------
 T mean solar angle (Greenwhich time) */
  days=jc2d*tj;
  astro_angles->sh_T=((days - (int)days)*24.0*15.0+180.0)*d2r;

/*------------------------------------------------------------------------------
 h mean solar Longitude */

  astro_angles->sh_h=(280.1895015 +omega_h *tj)*d2r;

/*------------------------------------------------------------------------------
 s mean lunar Longitude */

  astro_angles->sh_s=(277.0256206 +omega_s *tj)*d2r;

/*------------------------------------------------------------------------------
 p1 Longitude of solar perigee */

  astro_angles->sh_p1=(281.2208568 +omega_p1 *tj)*d2r;

/*------------------------------------------------------------------------------
 p Longitude of lunar perigee */

  astro_angles->sh_p=(334.3837214 +omega_p *tj)*d2r;

/*----------------------------------------------------------------------------*/
  astro_angles->sh_N =fmod(astro_angles->sh_N ,pi2);
  astro_angles->sh_s =fmod(astro_angles->sh_s ,pi2);
  astro_angles->sh_h =fmod(astro_angles->sh_h, pi2);
  radian_recale(&astro_angles->sh_p);
  astro_angles->sh_p1=fmod(astro_angles->sh_p1,pi2);

/*------------------------------------------------------------------------------
 I obliquity of the Moon's orbit (i.e. the angle with the equator)
 (Schureman, 1940, explanation of Tab6, p.156) */

  cosI=0.913694997 -0.035692561 *cos(astro_angles->sh_N);

  astro_angles->sh_I=acos(cosI);
  
  radian_recale(&astro_angles->sh_I);

/*------------------------------------------------------------------------------
 xi the longitude in the Moon's orbit, referenced at the intersection with the equator, of the vernal equinox
 (Schureman, 1940, explanation of Tab6) */

  sin2I=sin(astro_angles->sh_I);
  sh_tgn2=tan(astro_angles->sh_N/2.0);

  at1=atan(1.01883*sh_tgn2);
  at2=atan(0.64412*sh_tgn2);

  astro_angles->sh_xi=-at1-at2+astro_angles->sh_N;

  radian_recale(&astro_angles->sh_xi);

/*------------------------------------------------------------------------------
 nu the longitude of the intersection of the Moon's orbit and the Equator
 (Schureman, 1940, explanation of Tab6) */

  astro_angles->sh_nu=at1-at2;

/*------------------------------------------------------------------------------
 For constituents l2 k1 k2 */

  tgI2=tan(astro_angles->sh_I/2.0);
  
  /* Schureman (1940) eq.191 */
  P=astro_angles->sh_p-astro_angles->sh_xi;
  
  radian_recale(&P);
  
  /* Schureman (1940) eq.203 */
  q=atan2(0.483*sin(P),cos(P));
  /* Schureman (1940) eq.204 */
  astro_angles->sh_Qu=P-q;
  
  radian_recale(&astro_angles->sh_Qu);
  
#if 0
  FILE *f=fopen((string(__func__)+".dat").c_str(),"a");
  fprintf(f,"%g %g %g %g %g %g %g\n",astro_angles->t_julian,astro_angles->sh_p,astro_angles->sh_I,astro_angles->sh_xi,P,q,astro_angles->sh_Qu);
  fclose(f);
#endif
  
  /* Schureman (1940) eq.213 */
  t2=tgI2*tgI2;
  t4=t2*t2;
  astro_angles->sh_x1ra=sqrt(1.0-12.0*t2*cos(2.0*P)+36.0*t4);
  
  /* Schureman (1940) eq.214 */
  p=sin(2.0*P);
  q=1.0/(6.0*t2)-cos(2.0*P);
  astro_angles->sh_R=atan(p/q);
  
  /* Schureman (1940) eq.224 */
  p=sin(2.0*astro_angles->sh_I)*sin(astro_angles->sh_nu);
  q=sin(2.0*astro_angles->sh_I)*cos(astro_angles->sh_nu)+0.3347;
  astro_angles->sh_nuprim=atan(p/q);

  /* Schureman (1940) eq.232 */
  s2=sin(astro_angles->sh_I)*sin(astro_angles->sh_I);
  p=s2*sin(2.0*astro_angles->sh_nu);
  q=s2*cos(2.0*astro_angles->sh_nu)+0.0727;
  astro_angles->sh_nusec=0.5*atan(p/q);

/*----------------------------------------------------------------------------*/
  if (verbose) {
/*     STDOUT_BASE_LINE("%d/%d/%d \n",t_reference.day,t_reference.month,t_reference.year); */
    print_astro_angles(*astro_angles);
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_astro_angles(const astro_angles_t &astro_angles)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  STDOUT_BASE_LINE("s=%9.3f  h=%9.3f  p=%9.3f   p1=%9.3f \n",astro_angles.sh_s*r2d,astro_angles.sh_h*r2d,astro_angles.sh_p*r2d,astro_angles.sh_p1*r2d);
  STDOUT_BASE_LINE("I=%9.3f  N=%9.3f \n",astro_angles.sh_I*r2d,astro_angles.sh_N*r2d);
  STDOUT_BASE_LINE("xi=%9.3f nu=%9.3f nu'=%9.3f nu\"=%9.3f \n",astro_angles.sh_xi*r2d,astro_angles.sh_nu*r2d,astro_angles.sh_nuprim*r2d,astro_angles.sh_nusec*r2d);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void init_argument(astro_angles_t *astro_angles, date_t reference, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///wrapper for astronomic_angle(astro_angles_t *astro_angles, double t_julian, int verbose)
/**
\param reference model reference date and time
*/
/*----------------------------------------------------------------------------*/
{
  int N;
  double t_julian;                    //time elapsed between t_schureman and t_reference in julian centuries (36525 days)
  date_t t_schureman(1900,1,1,0.0);   //time reference for astronomical abaques
  spectrum_t spectrum;

  N=  julian_day(reference)
     -julian_day(t_schureman);

  t_julian=(N + reference.second/d2s) /jc2d;
  
  astronomic_angle(astro_angles,t_julian,verbose);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 void init_argument(astro_angles_t *astro_angles, double startd, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// wrapper for init_argument(astro_angles_t *astro_angles, date_t reference, int *verbose)
/**
\param startd seconds since 1 Jan 1950
*/
/*----------------------------------------------------------------------------*/
{
  date_t reference;
  ///calls poctime_getdatecnes()
  reference=poctime_getdatecnes(startd,'s');
  init_argument(astro_angles,reference,verbose);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 extern "C" void init_argument_(astro_angles_t *astro_angles, double *startd, int *verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Fortran wrapper for init_argument(astro_angles_t *,const double,int *)
/**
\param startd seconds since 1 Jan 1950
*/
/*----------------------------------------------------------------------------*/
{
  init_argument(astro_angles,*startd,*verbose);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double nodal_factor(const astro_angles_t &astro_angles, int formula)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///returns nodal amplitude factor
/**
\todo 2013-07-08 Florent Lyard and Damien Allain : non-linear waves's factor should be generated, removing a lot of code.
*/
/*----------------------------------------------------------------------------*/
{
  double s,f,P;
  double f1,f2;

  switch (formula) {

/* formule 0, solar waves */
    case 0:
      f=1.0;
      break;

/* formule 1, compound waves (M2 x M2, 78 x 78)*/
    case 1:
      f1=nodal_factor(astro_angles,78);
      f=f1*f1;
      break;

/* formule 2, compound waves (78 x 0) */
    case 2:
      f=nodal_factor(astro_angles,78);
      break;

/* formule 3, compound waves (M2 x K2, 78 x 235) */
    case 3:
      f1=nodal_factor(astro_angles,78);
      f2=nodal_factor(astro_angles,235);
      f=f1*f2;
      break;

/* formule 4, compound waves (M2 x M2 x M2, 78 x 78 x 78) */
    case 4:
      f1=nodal_factor(astro_angles,78);
      f=f1*f1*f1;
      break;

/* formule 5, compound waves (M2 x K1, 78 x 227) */
    case 5:
      f1=nodal_factor(astro_angles,78);
      f2=nodal_factor(astro_angles,227);
      f=f1*f2;
      break;

/* formule 6, compound waves (M2 x M2 x K1, 78 x 78 x 227) */
    case 6:
      f1=nodal_factor(astro_angles,78);
      f2=nodal_factor(astro_angles,227);
      f=f1*f1*f2;
      break;

/* formule 7, compound waves (O1 x ..., 75) */
    case 7:
      f=nodal_factor(astro_angles,75);
      break;

/* formule 8, compound waves (O1 x O1, 75 x 75)  */
    case 8:
      f1=nodal_factor(astro_angles,75);
      f2=nodal_factor(astro_angles,75);
      f=f1*f2;
      break;

/* formule 9, compound waves (M2 x M2 x K2, 78 x 78 x 235) */
    case 9:
      f1=nodal_factor(astro_angles,78);
      f2=nodal_factor(astro_angles,235);
      f=f1*f1*f2;
      break;
      
/* formule 10,  compound waves (78 x 227) */ // LR, add: 17/11/2009
    case 10:
      f = nodal_factor(astro_angles,78) * nodal_factor(astro_angles,227);
      break;
  
/* formule 11,  compound waves (75 x 0) */ // LR, add: 17/11/2009
    case 11:
      f = nodal_factor(astro_angles,75) * nodal_factor(astro_angles,0);
      break;
      
/* formule 12,  compound waves (78 x 78 x 78 x 0) */ // LR, add: 17/11/2009
    case 12:
      f = nodal_factor(astro_angles,78) * nodal_factor(astro_angles,78) * nodal_factor(astro_angles,78) * nodal_factor(astro_angles,0);
      break;
      
/* formule 13, compound waves (78 x 75)*/ // LR, add: 17/11/2009
    case 13:
      f = nodal_factor(astro_angles,78) * nodal_factor(astro_angles,75);
      break;
  
 /* formule 14, compound waves (235 x 0)  ===  (235)*/ // LR, add: 17/11/2009
    case 14:
      f = nodal_factor(astro_angles,235) * nodal_factor(astro_angles,0);
      break;
      
/* formule 15, compound waves (235 x 75) */ // LR, add: 17/11/2009
    case 15:
      f = nodal_factor(astro_angles,235) * nodal_factor(astro_angles,75);
      break;
      
/* formule 16, compound waves (78 x 0 x 0)  ===  (78)*/  // LR, add: 17/11/2009
    case 16:
      f = nodal_factor(astro_angles,78) * nodal_factor(astro_angles,0) * nodal_factor(astro_angles,0);
      break;
      
/* formule 17,  compound waves (227 x 0) */ // LR, add: 17/11/2009
    case 17:
      f = nodal_factor(astro_angles,227) * nodal_factor(astro_angles,0);
      break;
      
/* formule 18,  compound waves (78 x 78 x 78 ) */ // LR, add: 17/11/2009
    case 18:
      f = nodal_factor(astro_angles,78) * nodal_factor(astro_angles,78) * nodal_factor(astro_angles,78);
      break;
  
/* formule 19, compound waves (78 x 0 x 0 x 0)  ===  (78)*/ // LR, add: 17/11/2009
    case 19:
      f = nodal_factor(astro_angles,78) * nodal_factor(astro_angles,0) * nodal_factor(astro_angles,0) * nodal_factor(astro_angles,0);
      break;
      

/* All these are from Schureman (1940, from p25 onwards) #WAVELISTURL */

/* formule 73 */
    case 73:
      s=sin(astro_angles.sh_I);
      f=(2./3.-s*s)/0.5021;
      break;

/* formule 74 */
    case 74:
      s=sin(astro_angles.sh_I);
      f=s*s/0.1578;
      break;

/* formule 75 */
    case 75:
      s=cos (astro_angles.sh_I/2);
      f=sin (astro_angles.sh_I)*s*s/0.3800;
      break;

/* formule 76 */
    case 76:
      f=sin (2*astro_angles.sh_I)/0.7214;
      break;

/* formule 77 */
    case 77:
      s=sin (astro_angles.sh_I/2);
      f=sin (astro_angles.sh_I)*s*s/0.0164;
      break;

/* formule 78 */
    case 78:
      s=cos (astro_angles.sh_I/2);
      f=s*s*s*s/0.9154;
      break;

/* formule 79 */
    case 79:
      s=sin(astro_angles.sh_I);
      f=s*s/0.1565;
      break;

/* formule 141, p.36 */
    case 141:
      s = sin(astro_angles.sh_I);
      f=(s - 5./4. * s * s)/0.3192;
      break;

/* formule 144 */
    case 144:
      s = sin(astro_angles.sh_I/2);
      f=(1 - 10 * s * s + 15 * s * s * s * s) * square( cos(astro_angles.sh_I/2) )/0.5873;
      break;

/* formule 149 */
    case 149:
      s=cos (astro_angles.sh_I/2);
      f=s*s*s*s*s*s/0.8758;
      break;

/* formule 206, corrected */
    case 206:
      /* see Schureman (1940, eq.191 or eq.204) */
      P=astro_angles.sh_p-astro_angles.sh_xi;
      /* see Schureman (1940, eq.197) */
      s=sqrt(2.310+1.435*cos(2*P));
      /* see Schureman (1940, eq.207 and explanation of eq.206) */
      f1=nodal_factor(astro_angles,75);
      f=f1*s/1.52;
      break;

/* formule 215 */
    case 215:
      s=cos (astro_angles.sh_I/2);
      f=s*s*s*s/0.9154*astro_angles.sh_x1ra;
      break;

/* formule 227 */
    case 227:
      s=sin (2*astro_angles.sh_I);
      f=sqrt (0.8965*s*s+0.6001*s*cos (astro_angles.sh_nu)+0.1006);
      break;

/* formule 235 */
    case 235:
      s=sin (astro_angles.sh_I);
      f=sqrt (19.0444*s*s*s*s+2.7702*s*s*cos (2*astro_angles.sh_nu)+.0981);
      break;

    default:
      TRAP_ERR_EXIT(ENOEXEC,"formula %d not coded\n",formula);
    }

  return(f);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double nodal_phase(const astro_angles_t &astro_angles, tidal_wave w)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Returns the phase correction u due to nodal motion. Units are radian.

----------------------------------------------------------------------*/
{
  double u=0;

//bug 28-05-2010
//  u = astro_angles.sh_xi * w.nksi + astro_angles.sh_nu * w.nnu0 + astro_angles.sh_nuprim * w.nnu1 + astro_angles.sh_nusec * w.nnu2;
  u += astro_angles.sh_xi * w.nksi ;
  u += astro_angles.sh_nu * w.nnu0 ;
  u += astro_angles.sh_nuprim * w.nnu1 ;
  u += astro_angles.sh_nusec * w.nnu2 ;
  u += astro_angles.sh_R * w.Ra ;
  u += astro_angles.sh_Qu * w.Qu ;
  
  u = fmod(u,pi2);
  
  return(u);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double greenwhich_argument(const astro_angles_t &astro_angles, tidal_wave w)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** Gives the phase (in rad) of the tidal potential relative to the Greenwhich meridian
(e.g. the position of the fictuous celestial body).
\param w the wave
\returns the phase (in rad) of the tidal potential relative to the Greenwhich meridian
*/
/*----------------------------------------------------------------------------*/
{
  double V0;

  V0 = astro_angles.sh_T * w.nT + astro_angles.sh_s * w.ns + astro_angles.sh_h * w.nh + astro_angles.sh_p * w.np + astro_angles.sh_p1 * w.np1 + w.shift * d2r;
  V0 = fmod(V0,pi2);
  return(V0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double potential(const astro_angles_t &astro_angles, tidal_wave w, float x, float y)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Returns the phase of the tidal potential relative to the Greenwhich
  meridian (e.g. the position of the fictuous celestial body). Units are
  radian.

----------------------------------------------------------------------*/
{
  double P;

  P = w.Ap * cos(greenwhich_argument(astro_angles,w));

  return(P);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tidal_equilibrium(tidal_wave w, double x, double y, float *amp, float *pha)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

Compute astronomical potential constants. Parameters are :
  
  - w
  - x : longitude in degrees
  - y : latitude in degrees
  - *amp : will be set to amplitude in m
  - *pha : will be set to phase in degrees

Returns 0 if the potential is coded for and -1 otherwise.

------------------------------------------------------------------------------*/
{
  float a, G, yr, C, S;
  float k2 = 0.3, h2 = 0.6;

  a = w.Ap * (1 + k2 - h2);
  yr=y*d2r;
  
  switch(w.order){
  case 2:
  
    switch (w.nT) {
      case (0):

/*######################### Long period order 2 tide #########################

  From Schureman (1940) eq.62 p.21
   potential/g = A*(1/2 -3/2 sin^2(Y))*cos(w*t+V0)

----------------------------------------------------------------------*/

        C = cos(yr);
        S = sin(yr);
        a *= 0.5 - 1.5 * S * S;
        G = -0.;
        break;

      case (1):

/*########################### Diurnal order 2 tide ###########################

  From Schureman (1940) eq.63 p.21-22
   potential/g = A*sin(2Y)*cos(w*t+V0+X)

----------------------------------------------------------------------*/

        C = cos(yr);
        S = sin(yr);
        a *= 2 * S * C;
        G = -x;
        break;

      case (2):

/*######################### Semi-diurnal order 2 tide #########################

  From Schureman (1940) eq.64 p.22
   potential/g = A*cos^2(Y)*cos(w*t+V0+2*X)

----------------------------------------------------------------------*/

        C = cos(yr);
        S = sin(yr);
        a *= C * C;
        G = -2 * x;
        break;

/*####################### non-astronomical order 2 tide #######################
   potential/g = 0

----------------------------------------------------------------------*/
      default:
        *amp=-1.0;
        *pha= 0.0;
        return(-1);
        break;
      }
    
    break;
  
  case 3:
    switch (w.nT) {
      case (0):

/*######################### Long period order 3 tide #########################

  From Schureman (1940) eq.137 p.35
   potential/g = A*sin(Y)*(cos^2(Y)-2/5)*cos(w*t+V0)

----------------------------------------------------------------------*/

        C = cos(yr);
        S = sin(yr);
        a *= S*(C*C-2./5.);
        G = -0.;
        break;

      case (1):

/*########################### Diurnal order 3 tide ###########################

  From Schureman (1940) eq.138 p.35
   potential/g = A*cos(Y)*(cos^2(Y)-4/5)*cos(w*t+V0+X)

----------------------------------------------------------------------*/

        C = cos(yr);
        S = sin(yr);
        a *= C*(C*C-4./5.);
        G = -x;
        break;

      case (3):

/*######################### 3rd-diurnal order 3 tide #########################

  From Schureman (1940) eq.140 p.35
   potential/g = A*cos^3(Y)*cos(w*t+V0+2*X)

----------------------------------------------------------------------*/

        C = cos(yr);
        S = sin(yr);
        a *= C * C * C;
        G = -3 * x;
        break;

/*####################### non-astronomical tide #######################
   potential/g = 0

----------------------------------------------------------------------*/
      default:
        *amp=-1.0;
        *pha= 0.0;
        return(-1);
        break;
      }
    
    break;
  
  default:
    STDERR_BASE_LINE("Not coded for %dth order waves like %s\n",w.order,w.name);
    return -1;
    }

  //printf("tidal_equilibrium, %s deduced from static equilibrium\n",w.name);
  *amp = a;
  *pha = G;
  return(0);

}

typedef struct{const char *u;float f;} unit2factor_t;

const unit2factor_t
  ///time unit to hour factor table. \note Damien Allain 2012-05-07 : only used in omega_from_frequency_or_period()
  tu2ht[]={
    {"h",1.f},
    {"d",24.f},
    {"",NAN}},
   ///angle unit to degree factor table. \note Damien Allain 2012-05-07 : only used in omega_from_frequency_or_period()
  au2dt[]={
    {"deg",1.f},
    {"rad",r2d},
    /* no unit means cycle */
    {"",360.f}};


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float unit_to_factor(const char *unit,const unit2factor_t *table)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Finds the conversion factor of a unit in a given table
/**
\param *unit
\param *table
\returns the conversion factor
\note Damien Allain 2012-05-07 : only used in omega_from_frequency_or_period()
*/
/*----------------------------------------------------------------------------*/
{
  int i;
  
  for(i=0;i<3;i++){
    const unit2factor_t *ti=&table[i];
    if(ti->u[0]==0 || strncmp(table[i].u,unit)==0)
      return ti->f;
    }
  
  TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double omega_from_frequency_or_period(const float fOrT, const char *unit)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Gives a wave from its period or frequency
/**
\param fOrT period or frequency
\param *unit
\returns the wave

It uses unit_to_factor() with ::tu2ht and ::au2dt.
*/
/*----------------------------------------------------------------------------*/
{
  double omega=NAN,tf,af;
  char *unit_=strdup(unit);
  int dI;
  
/*---------------------------------------------------------------------*//**<h1>
  if \c unit is \c "dph", replace with \c "deg/h" </h1>*/
  if(strcmp(unit_,"dph")==0){
    free(unit_);
    unit_=strdup("deg/h");
    }
  
/*---------------------------------------------------------------------*//**<h1>
  detect from the presence of a \c '/' in \c unit whether it is a period or a frequency </h1>*/
  char *s=strchr(unit_,'/');
  
  if(s!=0){
    dI=s-unit_;
    //frequency
    char *angleUnit=strdup(unit_);
    const char *timeUnit=&unit_[dI+1];
    angleUnit[dI]='\0';
    
    af=unit_to_factor(angleUnit,au2dt);
    tf=unit_to_factor(timeUnit,tu2ht);
    omega=fOrT*af/tf;
    
    free(angleUnit);
    }
  else{
    //period
    tf=unit_to_factor(unit_,tu2ht);
    omega=360./(fOrT*tf);
    }
  
  free(unit_);
  
  return omega;
}


#define UNUSED_3RD_DOODSON_DIGIT "This does not take into account the 3rd from the end digit. Then again, no real wave seam to have it different from 5."


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  tidal_wave wave_from_doodson(int tau, int s, int h, int p, int N_, int p1, float coef, const char *name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

Gives a wave from its Doodson parameters.

N_ parameter is ignored.

See also wave_from_doodson(const int).

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  tidal_wave w=wNUL;
  
  if(name!=0)
    strcpy(w.name,name);
  
  w.nT=tau;
  w.ns=s-tau;
  w.nh=h+tau;
  ///\note 2011-10-20 Damien Allain : see #P_DISAGREEMENT
  w.np=p;
  ///\note 2011-10-11 Damien Allain : see #UNUSED_3RD_DOODSON_DIGIT
  //w.N'=N_;
  w.np1=p1;
  
  w.init();
  
  if(coef<=0.)
    w.shift=0;
  else
    w.shift=180;
  
  w.Ap=fabs(coef)/1.594;
  
  return w;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  tidal_wave wave_from_doodson(const int doodson)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Gives a wave from its IHO-like Doodson number (\b 7 digits for diurnal to 9th-diurnal)
/**
\param doodson Doodson number
\returns the wave

See also the IHO list of waves from #IHO_LIST_URL
*/
/*----------------------------------------------------------------------------*/
{
  tidal_wave w=wNUL;
  int r=doodson,f=1000000;//remainder,factor
  int tau, s, h, p, N_, p1;
  sprintf(w.name,"%d",doodson);
  
  tau=r/f;//7th from the end
  
  r%=f;f/=10;
  s=r/f-5;//6th from the end
  
  r%=f;f/=10;
  h=r/f-5;//5th from the end
  
  r%=f;f/=10;
  ///\note 2011-10-20 Damien Allain : see #P_DISAGREEMENT
  p=r/f-5;//4th from the end
  
  r%=f;f/=10;
  ///\note 2011-10-11 Damien Allain : see #UNUSED_3RD_DOODSON_DIGIT
  N_=r/f-5;//3rd from the end
  
  r%=f;f/=10;
  p1=r/f-5;//2nd from the end
  
  w=wave_from_doodson(tau, s, h, p, N_, p1, 0.f, 0);
  w.init();
  
  r%=f;f/=10;
  w.shift=(5-r/f)*90.;//1st from the end
  return w;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int doodson_from_wave(const tidal_wave & w)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Gives the IHO-like Doodson number (\b 7 digits for diurnal to 9th-diurnal) of a wave
/**
\param w the wave
\returns the Doodson number

See also the IHO list of waves from #IHO_LIST_URL
*/
/*----------------------------------------------------------------------------*/
{
  int r=0;//returned number
  r+=w.nT;//7th from the end
  r*=10;
  r+=w.ns+5+w.nT;//6th from the end
  r*=10;
  r+=w.nh+5-w.nT;//5th from the end
  r*=10;
  ///\note 2011-10-20 Damien Allain : see #P_DISAGREEMENT
  r+=w.np+5;//4th from the end
  r*=10;
  ///\note 2011-10-11 Damien Allain : see #UNUSED_3RD_DOODSON_DIGIT
  r+=5;//+w.N'//3rd from the end
  r*=10;
  r+=w.np1+5;//2nd from the end
  r*=10;
  r+=5-w.shift/90.;//1st from the end
  return r;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double tide_omega(tidal_wave wave, const char *units)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Gives the pulsation of the tide in degrees per hour
/**
\date reviewed 26 Jul 2011
\author Damien Allain

\param wave
\returns the pulsation of the tide in degrees per hour
*/
/*----------------------------------------------------------------------------*/
{
  double omega; //pulsation
  double scale;
  
  if(strcmp(units,"dph")==0) {
    scale = jc2d * 24.;               //conversion from degree per julian century to degree per hour
    }
  else if(strcmp(units,"rph")==0) {
    scale = jc2d * 24. *180.0 /M_PI;  //conversion from degree per julian century to radian per hour
    }
  else {
    scale = jc2d * 24.;               //conversion from degree per julian century to degree per hour
    }
  
  omega = omega_T * wave.nT + omega_s * wave.ns + omega_h * wave.nh +
          omega_p * wave.np + omega_p1 * wave.np1;
  omega /= scale;

  return (omega);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double tide_period(tidal_wave wave, bool initIf0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Gives the period of the tide in hours
/**
\date reviewed 26 Jul 2011
\author Damien Allain

\param wave
\returns the period of the tide in hours
*/
/*----------------------------------------------------------------------------*/
{
  double period;
  if(wave.omega==0. and initIf0){
    wave.init();
    }
  period=360./wave.omega;
  return (period);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void alias_frequency(double alias, double *omega)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// aliasing frequency
/**
\param alias aliasing frequency
\param[in,out] *omega frequency to be aliased, in the same frequency unit as \c f
*/
/*----------------------------------------------------------------------------*/
{
  const int
    n=round(*omega/alias);
  *omega=fabs(*omega-alias*n);
}


