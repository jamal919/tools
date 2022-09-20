#include "../config.h"

! MAKE SURE THE FILE BELOW IS UP TO DATE WITH libtools.a !
#include "tools-fortran-sizes.def"

#define STDOUT_FILE_LINE WRITE(0,*) __FILE__//":",__LINE__,":",


PROGRAM MAIN

use TTB

spectrum_t spectrum
!   tidal_wave *waves;     /*  array of tidal waves              */
!   int n;                /*  number of waves in array          */
!   int nmax;             /*  maximum number of waves in array  */

astro_angles_t angles

! TODO: SET THIS TO 1 IF YOUR MODEL DOES NODAL CORRECTIONS !!!
INTEGER, PARAMETER :: nodal=1 ! MARS does not do nodal corrections yet !

INTEGER nw ! number of waves

INTEGER nt !Number of time steps. Set this to 0 if you do not know them :)
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: time !time steps, if you know them.
REAL(KIND=8) t0,t,tL,dt !start time, current time, last time, time increment
! Wave names are considered to be 10 characters long.
CHARACTER(LEN=11), ALLOCATABLE, DIMENSION(:) :: waves !wave names
INTEGER ai ! atlas argument index
INTEGER wi ! wave argument index

INTEGER nij !buffer size
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: buffer
REAL(KIND=8), ALLOCATABLE, DIMENSION(:) :: lon,lat !MUST BE REAL(KIND=8)

!amplitudes and phases of constants
REAL(KIND=4), ALLOCATABLE, DIMENSION(:,:) :: amplitudes,degrees
REAL(KIND=4) mask
REAL(KIND=4) nan

INTEGER, PARAMETER :: outfile=100
INTEGER, PARAMETER :: infile=101

INTEGER narg ! number of arguments
CHARACTER(LEN=1024) :: arg0
CHARACTER(LEN=1024) :: controlPath,controlType,atlasPath
CHARACTER(LEN=1024) :: lonVar,latVar,ampVar,phaVar
INTEGER controlPathLen,controlPathPoint
INTEGER status,ncid,varid,ndims
INTEGER, ALLOCATABLE, DIMENSION(:) :: dims,dimlens
CHARACTER(LEN=80), ALLOCATABLE, DIMENSION(:) :: dimnames

nan=0
nan=nan/nan

narg=iargc()
WRITE(0,*) narg," arguments"
if (narg .lt. 3) then
  call getarg(0,arg0)
  WRITE(0,*) trim(arg0)," control.dat atlas.nc ampVar phaVar wave1 [wave2 ...]"
  WRITE(0,*) trim(arg0)," control.nc lonVar latVar atlas.nc ampVar phaVar wave1 [wave2 ...]"
  call exit(0)
  end if
 call getarg(1,controlPath)

! champs_Meteo.nc time boundaries
t0=1830297600 !2008/01/01 00:00:00 in CNES seconds
t0=1577836800 !2000/01/01 00:00:00 in CNES seconds
tL=1840665600 !2008/04/30 00:00:00 in CNES seconds
dt=3600
nt=(tL-t0)/dt+1
nt=8
WRITE(0,*) "number of time steps:",nt

! allocate the buffer to the number of control points
 controlPathLen=len(trim(controlPath))
 controlPathPoint=index(controlPath,".",1==1)+1
 controlType=controlPath(controlPathPoint:controlPathLen)
WRITE(0,*) "reading ",trim(controlType)," file ",trim(controlPath)

if (controlType=="dat") then
  OPEN(UNIT=infile,FILE=trim(controlPath),STATUS='OLD')
  READ(infile,*) nij
  allocate(lon(nij))
  allocate(lat(nij))
  do i=1,nij
    READ(infile,*) lon(i),lat(i)
    WRITE(0,*) i,lon(i),lat(i)
    end do
  CLOSE(infile)
  ai=2
  
else if (controlType=="nc") then
#if HAVE_LIBNETCDFF
  call getarg(2,lonVar)
  call getarg(3,latVar)
  ai=4
  status=nf_open(controlPath,0,ncid)
  call trap_error(status,'nf_open("'//trim(controlPath)//'",0,) error')
  status=nf_inq_varid(ncid,lonVar,varid)
  call trap_error(status,'nf_inq_varid(,"'//trim(lonVar)//'",) error')
  status=nf_inq_varndims(ncid,varid,ndims)
  call trap_error(status,'nf_inq_varndims() error')
  allocate(dims(ndims))
  allocate(dimnames(ndims))
  allocate(dimlens(ndims))
  status=nf_inq_vardimid(ncid,varid,dims)
  call trap_error(status,'nf_inq_vardimid() error')
  nij=1
  do i=1,ndims
    status=nf_inq_dim(ncid,dims(i),dimnames(i),dimlens(i))
    call trap_error(status,'nf_inq_dim() error')
    nij=nij*dimlens(i)
    end do
  WRITE(0,*) "reading ",nij," coordinates..."
  allocate(lon(nij))
  allocate(lat(nij))
  status=nf_get_var_double(ncid,varid,lon)
  call trap_error(status,'nf_get_var_double() error')
  call trap_error(status,'nf_get_var_double() error')
  status=nf_inq_varid(ncid,latVar,varid)
  call trap_error(status,'nf_inq_varid(,"'//trim(latVar)//'",) error')
  status=nf_get_var_double(ncid,varid,lat)
  call trap_error(status,'nf_get_var_double() error')
#else
  WRITE(*,*) "Unavailable type ",trim(controlType)
  call exit(8)
#endif
else 
  WRITE(*,*) "Unknown type ",trim(controlType)
  call exit(-1)
  end if

call getarg(ai,atlasPath)
call getarg(ai+1,ampVar)
call getarg(ai+2,phaVar)
wi=ai+2

! waves
! waves=(/("M2"),("S2"),("N2"),("K2"),("K1"),("O1"),("P1"),("Q1")/)
! waves=(/("Z0 "),("Q1 "),("O1 "),("P1 "),("K1 "),("N2 "),("M2 "),("S2 "),("K2 "),("L2 "),("M4 "),("MS4")/)
! waves=(/("M2 "),("S2 "),("K2 ")/)
! nw=size(waves) DOES NOT WORK WITH ifort
nw=narg-wi
!nw=1 !force single wave
WRITE(0,*) "number of waves:",nw
allocate(waves(nw))
do i=1,nw
  call getarg(i+wi,waves(i))
  end do

allocate(buffer(nij))
allocate(amplitudes(nij,nw))
allocate(degrees(nij,nw))

WRITE(0,*) "init_spectrum(,",nw,")"
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

WRITE(0,*) 'tide_atlas2positions("'//trim(atlasPath)//' ",...)'
mask=nf_fill_float
call tide_atlas2positions(atlasPath,spectrum,ampVar,phaVar,lon,lat,nij,nw,amplitudes,degrees,mask,1,status)
call trap_error(status,'tide_atlas2positions() error')
do wi=1,nw
  WRITE(*,*) waves(wi)
  if (nij<100) then
    WRITE(*,*) amplitudes(:,wi)
    WRITE(*,*) degrees(:,wi)
  else
    controlPath="interp-test.nc"
    WRITE(0,*) "saving amplitudes to:",trim(controlPath)
    status=nf_create(controlPath,nf_clobber,ncid)
    call trap_error(status,'nf_create("'//trim(controlPath)//'",0,) error')
    do i=1,ndims
      status=nf_def_dim(ncid,dimnames(i),dimlens(i),dims(i))
      call trap_error(status,'nf_def_dim() error')
      nij=nij*dimlens(i)
      end do
    status=nf_def_var(ncid,"Ha",nf_double,ndims,dims,varid)
    call trap_error(status,'nf_def_var(,"Ha",,,,) error')
!     status=nf_put_att_real(ncid,varid,"_FillValue",nf_double,1,mask)
!     call trap_error(status,'nf_put_att_real(,,"_FillValue",,,,) error')
    status=nf_enddef(ncid)
    call trap_error(status,'nf_enddef() error')
    status=nf_put_var_real(ncid,varid,amplitudes)
    call trap_error(status,'nf_put_var_real() error')
    call exit(0)
    end if
  end do

WRITE(0,*) "harmonic_prediction"
do i=1,nt
  t=dt*(i-1)
  call harmonic_prediction(angles,buffer,t,nij,spectrum,nw,amplitudes,degrees,nodal)
  if (nij<100) then
    WRITE(*,*) t/86400,buffer
  else
    WRITE(*,*) t/86400,buffer(j)
    end if
  end do

END
