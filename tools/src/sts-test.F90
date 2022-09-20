
! MAKE SURE THE FILE BELOW IS UP TO DATE WITH libtools.a !
#include "tools-fortran-sizes.def"

#define __OUT_FILE_LINE__ WRITE(0,*) __FILE__//":",__LINE__,":",


PROGRAM MAIN

use TTB

spectrum_t spectrum
!   tidal_wave *waves;     /*  array of tidal waves              */
!   int n;                /*  number of waves in array          */
!   int nmax;             /*  maximum number of waves in array  */

INTEGER nw ! number of waves

CHARACTER(LEN=11), ALLOCATABLE, DIMENSION(:) :: waves !wave names

! buffers for amplitudes and phases of elevation and U and V speeds
REAL(KIND=8), ALLOCATABLE, DIMENSION(:,:,:) :: e_a,e_g,u_a,u_g,v_a,v_g
INTEGER ni,nj !buffers sizes
INTEGER status

REAL(KIND=4) nan

nan=0
nan=nan/nan

! waves
! > grep -v ' ' data/OBC-FES2004/tides.obc
waves=(/("O1 "),("K1 "),("N2 "),("M2 "),("S2 ")/)

ni=422
nj=502

nw=size(waves)
!nw=1 !force single wave
allocate(e_a(ni,nj,nw))
allocate(e_g(ni,nj,nw))
allocate(u_a(ni,nj,nw))
allocate(u_g(ni,nj,nw))
allocate(v_a(ni,nj,nw))
allocate(v_g(ni,nj,nw))

WRITE(0,*) "init_spectrum"
call init_spectrum(spectrum,nw)
do i=1,nw
  WRITE(0,*) "add_wave_to_spectrum",i,"of",nw
  ! SPACE-TERMINATED STRING BECAUSE OF C TO FORTRAN MESS !
  ! You should do your own Fortran subroutine for this.
  ! Actually Fortran strings are space padded anyway :) so this works
  ! if the strings are enforced to a long enough name (like 11, see above):
  !call add_wave_to_spectrum(spectrum,waves(i))
  ! but adding an extra space will not hurt, so this is safer:
  call add_wave_to_spectrum(spectrum,waves(i)//" ")
  end do

WRITE(0,*) "spectral_tugo..."
call spectral_tugo("spectral-CQN1xCQP0.intg","","sts-fortran-test",spectrum, &
  ni,nj,nw,e_a,e_g,u_a,u_g,v_a,v_g,status)
WRITE(0,*) "spectral_tugo returned",status

END
