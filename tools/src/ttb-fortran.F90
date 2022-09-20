
! MAKE SURE THE FILE BELOW IS UP TO DATE WITH libtools.a !
#include "tools-fortran-sizes.def"

#define __OUT_FILE_LINE__ WRITE(0,*) __FILE__//":",__LINE__,":",


MODULE TTB
#include <netcdf.inc>
  
  PUBLIC harmonic_prediction,tide_atlas2positions,spectral_tugo
  
!******************************************************************************
INTERFACE harmonic_prediction
  
  SUBROUTINE harmonic_prediction4(astro_angles,buffer,time,n,s,nwaves,amplitudes,degrees,nodal)
    astro_angles_t,INTENT(in) :: astro_angles
    REAL(KIND=8),DIMENSION(n),INTENT(OUT) :: buffer
    REAL(KIND=8),INTENT(IN) :: time
    INTEGER,INTENT(IN) :: n,nwaves,nodal
    spectrum_t,INTENT(in) :: s
    REAL(KIND=4),DIMENSION(n,nwaves),INTENT(IN) :: amplitudes,degrees
  END SUBROUTINE
  
  SUBROUTINE harmonic_prediction8(astro_angles,buffer,time,n,s,nwaves,amplitudes,degrees,nodal)
    astro_angles_t,INTENT(in) :: astro_angles
    REAL(KIND=8),DIMENSION(n),INTENT(OUT) :: buffer
    REAL(KIND=8),DIMENSION(n),INTENT(IN) :: time
    INTEGER,INTENT(IN) :: n,nwaves,nodal
    spectrum_t,INTENT(in) :: s
    REAL(KIND=8),DIMENSION(n,nwaves),INTENT(IN) :: amplitudes,degrees
  END SUBROUTINE
  
END INTERFACE

!******************************************************************************
INTERFACE tide_atlas2positions
  MODULE PROCEDURE tide_atlas2positionsf4,tide_atlas2positionsf8
END INTERFACE

!******************************************************************************
INTERFACE spectral_tugo
  MODULE PROCEDURE spectral_tugof4,spectral_tugof8
END INTERFACE


!******************************************************************************
 CONTAINS
  
  !****************************************************************************
  SUBROUTINE spectral_tugof4(paramP,deltaP,outD,wl,ni,nj,nw,ha,hg,ua,ug,va,vg,status)
    
    CHARACTER(LEN=*),INTENT(in) :: paramP,deltaP,outD
    spectrum_t,INTENT(in) :: wl
    INTEGER,INTENT(IN) :: ni,nj,nw
    REAL(KIND=4),DIMENSION(ni,nj,nw),INTENT(OUT) :: ha,hg,ua,ug,va,vg
     INTEGER,INTENT(OUT) :: status
   
    ! GIVING LENGTH OF STRINGS BECAUSE FORTRAN STRINGS ARE SPACE TERMINATED !
    call spectral_tugo4(paramP,len_trim(paramP),deltaP,len_trim(deltaP), &
      outD,len_trim(outD),wl,ha,hg,ua,ug,va,vg,status)
    
  END SUBROUTINE
  
  !****************************************************************************
  SUBROUTINE spectral_tugof8(paramP,deltaP,outD,wl,ni,nj,nw,ha,hg,ua,ug,va,vg,status)
    
    CHARACTER(LEN=*),INTENT(in) :: paramP,deltaP,outD
    spectrum_t,INTENT(in) :: wl
    INTEGER,INTENT(IN) :: ni,nj,nw
    REAL(KIND=8),DIMENSION(ni,nj,nw),INTENT(OUT) :: ha,hg,ua,ug,va,vg
    INTEGER,INTENT(OUT) :: status
    
    ! GIVING LENGTH OF STRINGS BECAUSE FORTRAN STRINGS ARE SPACE TERMINATED !
    call spectral_tugo8(paramP,len_trim(paramP),deltaP,len_trim(deltaP), &
      outD,len_trim(outD),wl,ha,hg,ua,ug,va,vg,status)
    
  END SUBROUTINE
  
  
  !****************************************************************************
  SUBROUTINE trap_error(status,msg)
    
    INTEGER,INTENT(IN) :: status
    CHARACTER(LEN=*),INTENT(IN) :: msg
    
    if (status /= 0) then
      WRITE(0,*) msg,"(",status,trim(nf_strerror(status)),")"
      call exit(status)
      end if
    
  END SUBROUTINE
  
  
  !****************************************************************************
  SUBROUTINE tide_atlas2positionsf4(aPC,WaveList,av,gv,lon,lat,npositions,nwaves,a,G,mask,verbose,status)
    
    CHARACTER(LEN=*),INTENT(in) :: aPC
    spectrum_t,INTENT(in) :: WaveList
    CHARACTER(LEN=*),INTENT(in) :: av,gv
    INTEGER,INTENT(IN) :: npositions,nwaves
    REAL(KIND=8),DIMENSION(npositions),INTENT(IN) :: lon,lat
    
    REAL(KIND=4),DIMENSION(npositions,nwaves),INTENT(OUT) :: a,G
    REAL(KIND=4),INTENT(IN) :: mask
    
    INTEGER,INTENT(IN) :: verbose
    INTEGER,INTENT(OUT) :: status
    
    ! GIVING LENGTH OF STRINGS BECAUSE FORTRAN STRINGS ARE SPACE TERMINATED !
    call tide_atlas2positions4(aPC,len_trim(aPC),WaveList,av,len_trim(av), &
      gv,len_trim(gv),lon,lat,npositions,a,G,mask,verbose,status)
    
  END SUBROUTINE
  
  !****************************************************************************
  SUBROUTINE tide_atlas2positionsf8(aPC,WaveList,av,gv,lon,lat,npositions,nwaves,a,G,mask,verbose,status)
    
    CHARACTER(LEN=*),INTENT(in) :: aPC
    spectrum_t,INTENT(in) :: WaveList
    CHARACTER(LEN=*),INTENT(in) :: av,gv
    INTEGER,INTENT(IN) :: npositions,nwaves
    REAL(KIND=8),DIMENSION(npositions),INTENT(IN) :: lon,lat
    
    REAL(KIND=8),DIMENSION(npositions,nwaves),INTENT(OUT) :: a,G
    REAL(KIND=8),INTENT(IN) :: mask
    
    INTEGER,INTENT(IN) :: verbose
    INTEGER,INTENT(OUT) :: status
    
    ! GIVING LENGTH OF STRINGS BECAUSE FORTRAN STRINGS ARE SPACE TERMINATED !
    call tide_atlas2positions8(aPC,len_trim(aPC),WaveList,av,len_trim(av), &
      gv,len_trim(gv),lon,lat,npositions,a,G,mask,verbose,status)
    
  END SUBROUTINE

END MODULE
