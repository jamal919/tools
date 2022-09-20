
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
/*
MODULE eos
  !!======================================================================
  !!                     ***  MODULE  eos  ***
  !! All routines dealing with the Equation Of State of sea water
  !!=====================================================================
  !! History : 2.1  !  2004   : J.M. Molines : Original code ported
  !!                                           from NEMO
  !!           3.0    12/2010 : J.M. Molines : Doctor norm + Lic.
  !!----------------------------------------------------------------------
  !!----------------------------------------------------------------------
  !!   routines      : description
  !!   sigma0     : compute sigma-0 
  !!   eosbn2     : compute Brunt Vaissala Frequency
  !!   sigmai     : compute sigma-i ( refered to a depth given in argument
  !!   albet      : Compute the ratio alpha/beta ( Thermal/haline exapnsion)
  !!   beta       : compute beta (haline expension)
  !!----------------------------------------------------------------------
  IMPLICIT NONE

  PRIVATE
  PUBLIC :: sigma0
  PUBLIC :: eosbn2
  PUBLIC :: sigmai
  PUBLIC :: albet
  PUBLIC :: beta

  INTERFACE sigmai
     MODULE PROCEDURE sigmai_dep, sigmai_dep2d
  END INTERFACE


  !!----------------------------------------------------------------------
  !! CDFTOOLS_3.0 , MEOM 2011
  !! $Id: eos.f90 539 2011-07-11 10:33:35Z molines $
  !! Copyright (c) 2010, J.-M. Molines
  !! Software governed by the CeCILL licence (Licence/CDFTOOLSCeCILL.txt)
  !!----------------------------------------------------------------------

CONTAINS
*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double sigma0 ( float **ptem, float **psal, int kpi, int  kpj)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION sigma0  ***
    !!
    !! ** Purpose : Compute the in situ density (ratio rho/rau0) and the
    !!              potential volumic mass (Kg/m3) from potential temperature 
    !!              and salinity fields using an equation of state defined 
    !!              through the namelist parameter neos.
    !!
    !! ** Method  : Jackett and McDougall (1994) equation of state.
    !!              The in situ density is computed directly as a function of
    !!              potential temperature relative to the surface (the opa t
    !!              variable), salt and pressure (assuming no pressure variation
    !!              along geopotential surfaces, i.e. the pressure p in decibars
    !!              is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!              with pressure                      p        decibars
    !!                   potential temperature         t        deg celsius
    !!                   salinity                      s        psu
    !!                   reference volumic mass        rau0     kg/m**3
    !!                   in situ volumic mass          rho      kg/m**3
    !!                   in situ density anomalie      prd      no units
    !!
    !!----------------------------------------------------------------------
*/
//    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal  ! temperature and salinity
//    INTEGER(KIND=4),                  INTENT(in) :: kpi, kpj    ! dimension of 2D arrays
//    REAL(KIND=8), DIMENSION(kpi,kpj)             :: sigma0      ! returned value
{
    int ji, jj;
    double **zws;
    double zt, zs, zsr, zrau0=1000.;
    double zr1, zr2, zr3, zr4;
//    !!----------------------------------------------------------------------
//    zws = 0.0;
    double **sigma0;
    for( jj = 0; jj < kpj ; jj++) {
      for( ji = 0; ji < kpi ; ji++) {
         sigma0[jj][ji] = 0.;
         }
      }
    for( jj = 0; jj < kpj ; jj++) {
      for( ji = 0; ji < kpi ; ji++) {
         zws[jj][ji] = sqrt( fabs( psal[jj][ji] ) );
         }
      }

    for( jj = 0; jj < kpj ; jj++) {
//       !
      for( ji = 0; ji < kpi ; ji++) {

          zt  = ptem [jj][ji];        //    ! interpolated T
          zs  = psal [jj][ji];        //    ! interpolated S
          zsr = zws  [jj][ji];        //    ! square root of interpolated S

//          ! compute volumic mass pure water at atm pressure
          zr1 = ( ( ( ( 6.536332e-9*zt-1.120083e-6 )*zt+1.001685e-4)*zt 
               -9.095290e-3 )*zt+6.793952e-2 )*zt+999.842594;
//          ! seawater volumic mass atm pressure
          zr2= ( ( ( 5.3875e-9*zt-8.2467e-7 )*zt+7.6438e-5 ) *zt   
               -4.0899e-3 ) *zt+0.824493;
          zr3= ( -1.6546e-6*zt+1.0227e-4 ) *zt-5.72466e-3;
          zr4= 4.8314e-4;

//          ! potential volumic mass (reference to the surface)
          sigma0[jj][ji] = ( zr4*zs + zr3*zsr + zr2 ) *zs + zr1 - zrau0;
         }
      }

    return(0.0);
}
/*
  END FUNCTION sigma0
*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double sigmai_dep (float **ptem, float **psal, float **pref, int kpi, int  kpj)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*
    !! --------------------------------------------------------------------
    !! ** Purpose :   Compute the  density referenced to pref (ratio rho/rau0) 
    !!       from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter neos.
    !!
    !! ** Method  :
    !!       Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!              reference volumic mass        rau0     kg/m**3
    !!              in situ volumic mass          rho      kg/m**3
    !!              in situ density anomalie      prd      no units
    !! --------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi,kpj    ! dimension of 2D arrays
    REAL(KIND=4),                     INTENT(in) :: pref       ! reference pressure (meters or db)
    REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmai_dep             ! return value
*/
{
    double dpr4=4.8314e-4, dpd=-2.042967e-2 , dprau0 = 1000.e0;

    int ji, jj ;
    double   **dlrs;
    double dlt, dls   ;   
    double dla, dla1, dlaw, dlb, dlb1, dlbw, dlc, dle, dlk0, dlkw ;
    double dlrhop, dlr1, dlr2, dlr3, dlref;

//    dlref      = pref;
    double **sigmai_dep;

    for( jj = 0; jj < kpj ; jj++) {
      for( ji = 0; ji < kpi ; ji++) {
         dlrs[jj][ji] = sqrt( fabs( psal[jj][ji] ) );
         }
      }

    for( jj = 0; jj < kpj ; jj++) {
      for( ji = 0; ji < kpi ; ji++) {

//          ! Convert T and S to double precision.
          dlt = (double) (ptem[jj][ji]);
          dls = (double) (psal[jj][ji]);

//          ! Compute the volumic mass of pure water at atmospheric pressure.
          dlr1=((((6.536332e-9*dlt-1.120083e-6)
               *dlt+1.001685e-4)
               *dlt-9.095290e-3)
               *dlt+6.793952e-2)
               *dlt+999.842594;

//          ! Compute the seawater volumic mass at atmospheric pressure.
          dlr2=(((5.3875e-9*dlt-8.2467e-7)
               *dlt+7.6438e-5)
               *dlt-4.0899e-3)
               *dlt+0.824493;

          dlr3=(-1.6546e-6*dlt+1.0227e-4)
               *dlt-5.72466e-3;

//          ! Compute the potential volumic mass (referenced to the surface).
          dlrhop=(dpr4*dls+dlr3*dlrs[jj][ji]+dlr2)*dls+dlr1;

//          ! Compute the compression terms.
          dle=(-3.508914e-8*dlt-1.248266e-8)
               *dlt-2.595994e-6;

          dlbw=(1.296821e-6*dlt-5.782165e-9)
               *dlt+1.045941e-4;

          dlb=dlbw+dle*dls;

          dlc=(-7.267926e-5*dlt+2.598241e-3 )
               *dlt+0.1571896;

          dlaw=((5.939910e-6*dlt+2.512549e-3)
               *dlt-0.1028859)
               *dlt-4.721788;

          dla=(dpd*dlrs[jj][ji]+dlc)*dls+dlaw;

          dlb1=(-0.1909078*dlt+7.390729)
               *dlt-55.87545;

          dla1=((2.326469e-3*dlt+1.553190)
               *dlt-65.00517)
               *dlt+1044.077;

          dlkw=(((-1.361629e-4*dlt-1.852732e-2)
               *dlt-30.41638)
               *dlt+2098.925)
               *dlt+190925.6;

          dlk0=(dlb1*dlrs[jj][ji]+dla1)*dls+dlkw;

//          ! Compute the potential density anomaly.
          sigmai_dep[jj][ji]=dlrhop/(1.0-dlref/(dlk0-dlref*(dla-dlref*dlb)))
               -dprau0;

        }
      }

    return(0.0);
}
/*
  END FUNCTION sigmai_dep

  FUNCTION sigmai_dep2d ( ptem, psal, pref, kpi,kpj)
    !! --------------------------------------------------------------------
    !! ** Purpose :   Compute the  density referenced to pref (ratio rho/rau0) 
    !!       from potential temperature and
    !!      salinity fields using an equation of state defined through the
    !!     namelist parameter neos.
    !!
    !! ** Method  :
    !!       Jackett and McDougall (1994) equation of state.
    !!         the in situ density is computed directly as a function of
    !!         potential temperature relative to the surface (the opa t
    !!         variable), salt and pressure (assuming no pressure variation
    !!         along geopotential surfaces, i.e. the pressure p in decibars
    !!         is approximated by the depth in meters.
    !!              prd(t,s,p) = ( rho(t,s,p) - rau0 ) / rau0
    !!              rhop(t,s)  = rho(t,s,0)
    !!         with pressure                      p        decibars
    !!              potential temperature         t        deg celsius
    !!              salinity                      s        psu
    !!              reference volumic mass        rau0     kg/m**3
    !!              in situ volumic mass          rho      kg/m**3
    !!              in situ density anomalie      prd      no units
    !! --------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    INTEGER(KIND=4),                  INTENT(in) :: kpi,kpj    ! dimension of 2D arrays
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: pref       ! reference pressure (meters or db) (2d Array)
    REAL(KIND=8), DIMENSION(kpi,kpj) :: sigmai_dep2d           ! return value

    REAL(kind=8), PARAMETER :: dpr4=4.8314d-4, dpd=-2.042967d-2 , dprau0 = 1000.d0

    INTEGER(KIND=4) :: ji, jj 
    REAL(KIND=8), DIMENSION (kpi,kpj) ::  dlrs
    REAL(KIND=8) :: dlt, dls      
    REAL(KIND=8) :: dla, dla1, dlaw, dlb, dlb1, dlbw, dlc, dle, dlk0, dlkw 
    REAL(kind=8) :: dlrhop, dlr1, dlr2, dlr3, dlref

    sigmai_dep2d = 0.d0
    DO jj = 1, kpj
       DO ji = 1, kpi
          dlrs(ji,jj) = SQRT( ABS( psal(ji,jj) ) )
       END DO
    END DO

    DO jj=1,kpj
       DO ji=1,kpi

          ! Convert T and S to double precision.
          dlt   = DBLE(ptem(ji,jj))
          dls   = DBLE(psal(ji,jj))
          dlref = DBLE(pref(ji,jj))

          ! Compute the volumic mass of pure water at atmospheric pressure.
          dlr1=((((6.536332d-9*dlt-1.120083d-6)&
               *dlt+1.001685d-4)&
               *dlt-9.095290d-3)&
               *dlt+6.793952d-2)&
               *dlt+999.842594d0

          ! Compute the seawater volumic mass at atmospheric pressure.
          dlr2=(((5.3875d-9*dlt-8.2467d-7)&
               *dlt+7.6438d-5)&
               *dlt-4.0899d-3)&
               *dlt+0.824493d0

          dlr3=(-1.6546d-6*dlt+1.0227d-4)&
               *dlt-5.72466d-3

          ! Compute the potential volumic mass (referenced to the surface).
          dlrhop=(dpr4*dls+dlr3*dlrs(ji,jj)+dlr2)*dls+dlr1

          ! Compute the compression terms.
          dle=(-3.508914d-8*dlt-1.248266d-8)&
               *dlt-2.595994d-6

          dlbw=(1.296821d-6*dlt-5.782165d-9)&
               *dlt+1.045941d-4

          dlb=dlbw+dle*dls

          dlc=(-7.267926d-5*dlt+2.598241d-3 )&
               *dlt+0.1571896d0

          dlaw=((5.939910d-6*dlt+2.512549d-3)&
               *dlt-0.1028859d0)&
               *dlt-4.721788d0

          dla=(dpd*dlrs(ji,jj)+dlc)*dls+dlaw

          dlb1=(-0.1909078d0*dlt+7.390729d0)&
               *dlt-55.87545d0

          dla1=((2.326469d-3*dlt+1.553190d0)&
               *dlt-65.00517d0)&
               *dlt+1044.077d0

          dlkw=(((-1.361629d-4*dlt-1.852732d-2)&
               *dlt-30.41638d0)&
               *dlt+2098.925d0)&
               *dlt+190925.6d0

          dlk0=(dlb1*dlrs(ji,jj)+dla1)*dls+dlkw

          ! Compute the potential density anomaly.
          sigmai_dep2d(ji,jj)=dlrhop/(1.0d0-dlref/(dlk0-dlref*(dla-dlref*dlb)))&
               -dprau0

       ENDDO
    ENDDO

  END FUNCTION sigmai_dep2d


  FUNCTION eosbn2 ( ptem, psal, pdep, pe3w, kpi, kpj, kup, kdown)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION eosbn2  ***
    !!
    !! ** Purpose :  Compute the local Brunt-Vaisala frequency at the time-
    !!               step of the input arguments
    !!
    !! ** Method :  UNESCO sea water properties
    !!             The brunt-vaisala frequency is computed using the
    !!              polynomial expression of McDougall (1987):
    !!              N^2 = grav * beta * ( alpha/beta*dk[ t ] - dk[ s ] )/e3w
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj,2), INTENT(in) :: ptem, psal ! temperaature salinity
    REAL(KIND=4)                                   :: pdep       ! reference depth
    REAL(KIND=4), DIMENSION(kpi,kpj),   INTENT(in) :: pe3w       ! e3w of the current layer
    INTEGER(KIND=4),                    INTENT(in) :: kpi, kpj   ! size of the array
    INTEGER(KIND=4),                    INTENT(in) :: kup, kdown ! index of levels up and down
    REAL(KIND=4), DIMENSION(kpi,kpj)               :: eosbn2     ! returned values

    INTEGER(KIND=4) :: ji, jj         ! dummy loop indices
    REAL(KIND=8)    :: zgde3w, zt, zs, zh
    REAL(KIND=8)    :: zalbet, zbeta
    REAL(KIND=8)    :: zgrav=9.81
    !!----------------------------------------------------------------------

    zh = pdep
    DO jj = 1, kpj
       DO ji = 1, kpi
          zgde3w = zgrav / pe3w(ji,jj)
          zt = 0.5 * ( ptem(ji,jj,kup) + ptem(ji,jj,kdown) )          ! potential temperature at w-point
          zs = 0.5 * ( psal(ji,jj,kup) + psal(ji,jj,kdown) ) - 35.0   ! salinity anomaly (s-35) at w-point

          zalbet = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt     &   ! ratio alpha/beta
               &                               - 0.203814e-03 ) * zt   &
               &                               + 0.170907e-01 ) * zt   &
               &   + 0.665157e-01                                      &
               &   +     ( - 0.678662e-05 * zs                         &
               &           - 0.846960e-04 * zt + 0.378110e-02 ) * zs   &
               &   +   ( ( - 0.302285e-13 * zh                         &
               &           - 0.251520e-11 * zs                         &
               &           + 0.512857e-12 * zt * zt           ) * zh   &
               &           - 0.164759e-06 * zs                         &
               &        +(   0.791325e-08 * zt - 0.933746e-06 ) * zt   &
               &                               + 0.380374e-04 ) * zh

          zbeta  = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt        &   ! beta
               &                            - 0.301985e-05 ) * zt      &
               &   + 0.785567e-03                                      &
               &   + (     0.515032e-08 * zs                           &
               &         + 0.788212e-08 * zt - 0.356603e-06 ) * zs     &
               &   +(  (   0.121551e-17 * zh                           &
               &         - 0.602281e-15 * zs                           &
               &         - 0.175379e-14 * zt + 0.176621e-12 ) * zh     &
               &                             + 0.408195e-10   * zs     &
               &     + ( - 0.213127e-11 * zt + 0.192867e-09 ) * zt     &
               &                             - 0.121555e-07 ) * zh

          eosbn2(ji,jj) = zgde3w * zbeta                                         &   ! N^2
               &          * ( zalbet * ( ptem(ji,jj,kup) - ptem(ji,jj,kdown) )   &
               &                     - ( psal(ji,jj,kup) - psal(ji,jj,kdown) ) )
       END DO
    END DO

  END FUNCTION eosbn2


  FUNCTION albet(  ptem, psal, pdep, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION  albet  ***
    !!
    !! ** Purpose :  Compute the ratio alpha/beta 
    !!
    !! ** Method  :  Follow Mc Dougal et al as in other functions
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj), INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4),                     INTENT(in) :: pdep       ! refererence depth
    INTEGER(KIND=4),                  INTENT(in) :: kpi, kpj   ! size of the arrays

    REAL(KIND=8), DIMENSION(kpi,kpj)             :: albet      ! returned value

    INTEGER(KIND=4) :: ji, jj      ! dummy loop index
    REAL(KIND=8)    :: zt, zs, zh  ! working local variables
    !!----------------------------------------------------------------------
    zh = pdep
    DO ji=1,kpi
       DO jj=1,kpj
          zt =  ptem(ji,jj)         ! potential temperature
          zs =  psal(ji,jj)- 35.0   ! salinity anomaly (s-35)

          albet(ji,jj) = ( ( ( - 0.255019e-07 * zt + 0.298357e-05 ) * zt   &   ! ratio alpha/beta
               &                               - 0.203814e-03 ) * zt   &
               &                               + 0.170907e-01 ) * zt   &
               &   + 0.665157e-01                                      &
               &   +     ( - 0.678662e-05 * zs                         &
               &           - 0.846960e-04 * zt + 0.378110e-02 ) * zs   &
               &   +   ( ( - 0.302285e-13 * zh                         &
               &           - 0.251520e-11 * zs                         &
               &           + 0.512857e-12 * zt * zt           ) * zh   &
               &           - 0.164759e-06 * zs                         &
               &        +(   0.791325e-08 * zt - 0.933746e-06 ) * zt   &
               &                               + 0.380374e-04 ) * zh
       END DO
    END DO

  END FUNCTION albet


  FUNCTION beta (  ptem, psal, pdep, kpi, kpj)
    !!---------------------------------------------------------------------
    !!                  ***  FUNCTION beta  ***
    !!
    !! ** Purpose :  Compute the beta 
    !!
    !! ** Method  :  Follow Mc Dougal et al as in other functions 
    !!
    !!----------------------------------------------------------------------
    REAL(KIND=4), DIMENSION(kpi,kpj),INTENT(in) :: ptem, psal ! temperature salinity
    REAL(KIND=4),                    INTENT(in) :: pdep       ! reference depth
    INTEGER(KIND=4),                 INTENT(in) :: kpi, kpj   ! size of the array
    REAL(KIND=8), DIMENSION(kpi,kpj)            :: beta       ! returned values

    INTEGER(KIND=4) :: ji, jj      ! dummy loop index
    REAL(KIND=8)    :: zt, zs, zh  ! working variables
    !!----------------------------------------------------------------------
    zh = pdep
    DO ji=1,kpi
       DO jj=1,kpj
          zt =  ptem(ji,jj)         ! potential temperature
          zs =  psal(ji,jj)- 35.0   ! salinity anomaly (s-35)

          beta(ji,jj)  = ( ( -0.415613e-09 * zt + 0.555579e-07 ) * zt      &   ! beta
               &                            - 0.301985e-05 ) * zt      &
               &   + 0.785567e-03                                      &
               &   + (     0.515032e-08 * zs                           &
               &         + 0.788212e-08 * zt - 0.356603e-06 ) * zs     &
               &   +(  (   0.121551e-17 * zh                           &
               &         - 0.602281e-15 * zs                           &
               &         - 0.175379e-14 * zt + 0.176621e-12 ) * zh     &
               &                             + 0.408195e-10   * zs     &
               &     + ( - 0.213127e-11 * zt + 0.192867e-09 ) * zt     &
               &                             - 0.121555e-07 ) * zh
       END DO
    END DO

  END FUNCTION beta

END MODULE eos
*/