PROGRAM MAIN

#define __OUT_FILE_LINE__ WRITE(0,*) __FILE__//":",__LINE__,":",

! libtools.a is made to depend to the automatically generated file below
#include "tools-fortran-sizes.def"

spectrum_t spectrum
!   tidal_wave *waves;     /*  array of tidal waves              */
!   int n;                /*  number of waves in array          */
!   int nmax;             /*  maximum number of waves in array  */

harmonic_t harmonic
!   double *A;///< harmonic matrix : [neq*neq]
!   double ***rhs;///< right-hand side (RHS) vector : [nrhs][nndes][neq]
!   double **cs, **sn; ///< cosines and sines : [spectrum.n][nframe]
!   int *pivot;
!   int nframes;///< number of frames
!   int neq;///< =2*spectrum.n
!   int nrhs;///< number of variables ALWAYS 1 FOR YOU
!   int nndes;///< number of points in the grid
!   spectrum_t spectrum;
!   bool *mask;

astro_angles_t angles

! TODO: SET THIS TO 1 IF YOUR MODEL DOES NODAL CORRECTIONS !!!
INTEGER, PARAMETER :: nodal=1 ! MARS does not do nodal corrections yet !

INTEGER nw ! number of waves
INTEGER ni,nj !grid size
INTEGER nt !Number of time steps. Set this to 0 if you do not know them :)
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: time !time steps, if you know them.
REAL(KIND=8) t0,t,tL,dt !start time, current time, last time, time increment
! Wave names are considered to be 10 characters long.
CHARACTER(LEN=11), ALLOCATABLE, DIMENSION(:) :: waves !wave names
INTEGER nij,neq !buffer size:=ni*nj. neq=2*nw.

REAL(KIND=4), ALLOCATABLE, DIMENSION(:) :: buffer
CHARACTER(LEN=10) :: nijs !format string
INTEGER dosave,hassaved
!harmonic matrix and right-hand-side vector
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: A,rhs,readA,readrhs
!real and imaginary parts of constants
REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: rs,is

REAL(KIND=4), ALLOCATABLE, DIMENSION(:) :: readbuf
INTEGER, PARAMETER :: outfile=100
CHARACTER(LEN=100) A_and_rhs_path
INTEGER, PARAMETER :: infile=101

! champs_Meteo.nc dimensions, for example :
ni = 1190
nj = 1392

! champs_Meteo.nc time boundaries
t0=1830297600 !2008/01/01 00:00:00 in CNES seconds
tL=1840665600 !2008/04/30 00:00:00 in CNES seconds
dt=3600
nt=(tL-t0)/dt+1 ! Set this to 0 if
nt=0   ! you do not know them.
WRITE(0,*) "number of known time steps:",nt

! waves
! waves=(/("M2"),("S2"),("N2"),("K2"),("K1"),("O1"),("P1"),("Q1")/)
waves=(/("Z0 "),("Q1 "),("O1 "),("P1 "),("K1 "),("N2 "),("M2 "),("S2 "),("K2 "),("L2 "),("M4 "),("MS4")/)

! allocate the buffers
nij=ni*nj ! to the number of model points
! or to the number of control points
OPEN(UNIT=infile,FILE='control.dat',STATUS='OLD')
READ(infile,*) nij
CLOSE(infile)

allocate(buffer(nij))
nw=size(waves)
!nw=1 !force single wave
allocate(rs(nij,nw))
allocate(is(nij,nw))
neq=nw*2
allocate(A(neq*neq))
allocate(rhs(nij*neq))

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

WRITE(0,*) "init_argument"
call init_argument(angles,t0,0)
WRITE(0,*) "harmonic_start",nij
call harmonic_start(harmonic,spectrum,nij,nt)

OPEN(UNIT=infile,FILE='series.dat',STATUS='OLD')
allocate(readbuf(nij*3))

hassaved=0
if(nt>0)then
  ! you know your times already
  allocate(time(nt))
  do i=1,nt
    time(i)=t0+dt*(i-1)
    enddo
  WRITE(0,*) "harmonic_init"
  ! give CNES time
  call harmonic_init(harmonic,time,nodal,angles)
  WRITE(0,*) "harmonic_storage?_i"
  !DO NOT OMP PARALLEL THIS LOOP
  do i=1,nt
    ! vvvvvvvvvvvv DO YOUR MODELLING BELOW vvvvvvvvvvvvvvvv
    READ(infile,*) (readbuf(j),j=1,nij*3)
    t=readbuf(1)*86400.+t0
    do j=1,nij
      buffer(j)=readbuf((j-1)*3+2)
      end do
    ! ...
    ! ^^^^^^^^^^^^ DO YOUR MODELLING ABOVE ^^^^^^^^^^^^^^^^
    call harmonic_storage4_i(harmonic,i,buffer)
    ! vvvvvvvvvvvv DO YOUR MODELLING BELOW vvvvvvvvvvvvvvvv
    ! ...
    ! ^^^^^^^^^^^^ DO YOUR MODELLING ABOVE ^^^^^^^^^^^^^^^^
    end do
else
  WRITE(0,*) "harmonic_storage?_t"
  !DO NOT OMP PARALLEL THIS LOOP
  do while (t<tL-dt*.5)
    ! vvvvvvvvvvvv DO YOUR MODELLING BELOW vvvvvvvvvvvvvvvv
    READ(infile,*) (readbuf(j),j=1,nij*3)
    t=readbuf(1)*86400.+t0
    do j=1,nij
      buffer(j)=readbuf((j-1)*3+2)
      end do
    ! ...
    dosave=1 ! set this to 1 to test post-processing
    ! ...
    ! ^^^^^^^^^^^^ DO YOUR MODELLING ABOVE ^^^^^^^^^^^^^^^^
    ! give time since time given to init_argument
    call harmonic_storage4_t(harmonic,t-t0,buffer,nodal,angles)
    if(dosave>0)then
      call harmonic_get_A_and_rhs(harmonic,A,rhs)
      hassaved=hassaved+1
      WRITE(A_and_rhs_path,*) hassaved
      OPEN(UNIT=outfile,FILE=trim(adjustl(A_and_rhs_path)),STATUS='REPLACE')
      WRITE(outfile,*) A,rhs
      CLOSE(UNIT=outfile)
      !zero the harmonic matrix
      A(:)=0
      rhs(:)=0
      call harmonic_set_A_and_rhs(harmonic,A,rhs)
      dosave=0
      end if
    ! vvvvvvvvvvvv DO YOUR MODELLING BELOW vvvvvvvvvvvvvvvv
    ! ...
    ! t=...
    ! ...
    ! ^^^^^^^^^^^^ DO YOUR MODELLING ABOVE ^^^^^^^^^^^^^^^^
    end do
  end if

CLOSE(UNIT=infile)

if(hassaved>0)then
  WRITE(*,*) "post-processing"
  allocate(readA(neq*neq))
  allocate(readrhs(nij*neq))
  A(:)=0
  rhs(:)=0
  do i=1,hassaved
    ! read your files
    WRITE(A_and_rhs_path,*) i
    OPEN(UNIT=outfile,FILE=trim(adjustl(A_and_rhs_path)),STATUS='OLD')
    READ(outfile,*) readA,readrhs
    CLOSE(UNIT=outfile)
    ! sum your files
    A(:)=A(:)+readA(:)
    rhs(:)=rhs(:)+readrhs(:)
    end do
  call harmonic_set_A_and_rhs(harmonic,A,rhs)
  endif

WRITE(*,*) "harmonic_analysis"
call harmonic_analysis(harmonic,tL-t0,1,rs,is)
WRITE(*,*) "real parts ###################################################"
WRITE(*,*) rs
WRITE(*,*) "imag parts ###################################################"
WRITE(*,*) is
WRITE(*,*) "moduluses ####################################################"
WRITE(*,*) sqrt(rs*rs+is*is)

END
