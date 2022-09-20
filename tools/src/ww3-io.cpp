
/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#include <config.h>

#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <stdarg.h>
#include <cmath>
#include <map>

#include <iostream>
#include <string>
#include <fstream>
#include <sstream>
#include <vector>
#include <iterator>

#include "tools-structures.h"
#include "map.h"

#include "map.h"

/* *----------------------------------------------------------------------------
$ Grid name (C*30, in quotes)
$ 'IROISE'

$
$ Frequency increment factor and first frequency (Hz) ---------------- $
$ number of frequencies (wavenumbers) and directions, relative offset
$ of first direction in terms of the directional increment [-0.5,0.5].
$ In versions 1.18 and 2.22 of the model this value was by definiton 0,
$ it is added to mitigate the GSE for a first order scheme. Note that
$ this factor is IGNORED in the print plots in ww3_outp.
$
$  1.1  0.04118  25  24  0.

$ Set model flags
$  - FLDRY         Dry run (input/output only, no calculation).
$  - FLCX, FLCY    Activate X and Y component of propagation.
$  - FLCTH, FLCK   Activate direction and wavenumber shifts.
$  - FLSOU         Activate source terms.
$
$   F T T T F T

$ Set time steps ----------------------------------------------------- $
$ - Time step information (this information is always read)
$     maximum global time step, maximum CFL time step for x-y and
$     k-theta, minimum source term time step (all in seconds).
$   60. 95. 60. 20.

$ Start of namelist input section ------------------------------------ $
$   Starting with WAVEWATCH III version 2.00, the tunable parameters
$   for source terms, propagation schemes, and numerics are read using
$   namelists. Any namelist found in the folowing sections up to the
$   end-of-section identifier string (see below) is temporarily written
$   to ww3_grid.scratch, and read from there if necessary. Namelists
$   not needed for the given switch settings will be skipped
$   automatically, and the order of the namelists is immaterial.
$   As an example, namelist input to change SWELLF and ZWND in the
$   Tolman and Chalikov input would be
$
$   &SIN2 SWELLF = 0.1, ZWND = 15. /
$
$ Define constants in source terms ----------------------------------- $
$
$ Stresses - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$   TC 1996 with cap    : Namelist FLX3
$                           CDMAX  : Maximum allowed CD (cap)
$                           CTYPE  : Cap type :
$                                     0: Discontinuous (default).
$                                     1: Hyperbolic tangent.
$
$ Linear input - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$   Cavaleri and M-R    : Namelist SLN1
$                           CLIN   : Proportionality constant.
$                           RFPM   : Factor for fPM in filter.
$                           RFHF   : Factor for fh in filter.
$
$ Exponential input  - - - - - - - - - - - - - - - - - - - - - - - - -
$   WAM-3               : Namelist SIN1
$                           CINP   : Proportionality constant.
$
$   Tolman and Chalikov : Namelist SIN2
$                           ZWND   : Height of wind (m).
$                           SWELLF : swell factor in (n.nn).
$                           STABSH, STABOF, CNEG, CPOS, FNEG :
$                                    c0, ST0, c1, c2 and f1 in . (n.nn)
$                                    through (2.65) for definition of
$                                    effective wind speed (!/STAB2).
$   WAM4 and variants  : Namelist SIN3
$                           ZWND    : Height of wind (m).
$                           ALPHA0  : minimum value of Charnock coefficient
$                           Z0MAX   : maximum value of air-side roughness z0
$                           BETAMAX : maximum value of wind-wave coupling
$                           SINTHP  : power of cosine in wind input
$                           ZALP    : wave age shift to account for gustiness
$                       TAUWSHELTER : sheltering of short waves to reduce u_star
$                         SWELLFPAR : choice of swell attenuation formulation
$                                    (1: TC 1996, 3: ACC 2008)
$                            SWELLF : swell attenuation factor
$     Extra parameters for SWELLFPAR=3 only
$                  SWELLF2, SWELLF3 : swell attenuation factors
$                           SWELLF4 : Threshold Reynolds number for ACC2008
$                  	    SWELLF5 : Relative viscous decay below threshold
$                             Z0RAT : roughness for oscil. flow / mean flow
$
$ Nonlinear interactions - - - - - - - - - - - - - - - - - - - - - - -
$   Discrete I.A.       : Namelist SNL1
$                           LAMBDA : Lambda in source term.
$                           NLPROP : C in sourc term. NOTE : default
$                                    value depends on other source
$                                    terms selected.
$                           KDCONV : Factor before kd in Eq. (n.nn).
$                           KDMIN, SNLCS1, SNLCS2, SNLCS3 :
$                                    Minimum kd, and constants c1-3
$                                    in depth scaling function.
$   Exact interactions  : Namelist SNL2
$                           IQTYPE : Type of depth treatment
$                                     1 : Deep water
$                                     2 : Deep water / WAM scaling
$                                     3 : Shallow water
$                           TAILNL : Parametric tail power.
$                           NDEPTH : Number of depths in for which
$                                    integration space is established.
$                                    Used for IQTYPE = 3 only
$                         Namelist ANL2
$                           DEPTHS : Array with depths for NDEPTH = 3
$
$ Dissipation  - - - - - - - - - - - - - - - - - - - - - - - - - - - -
$   WAM-3               : Namelist SDS1
$                           CDIS, APM : As in source term.
$
$   Tolman and Chalikov : Namelist SDS2
$                           SDSA0, SDSA1, SDSA2, SDSB0, SDSB1, PHIMIN :
$                                    Constants a0, a1, a2, b0, b1 and
$                                    PHImin.
$
$   WAM4 and variants   : Namelist SDS3
$                           SDSC1    : WAM4 Cds coeffient
$                           MNMEANP, WNMEANPTAIL : power of wavenumber
$                                    for mean definitions in Sds and tail
$                           SDSDELTA1, SDSDELTA2 : relative weights
$                                    of k and k^2 parts of WAM4 dissipation
$                           SDSLF, SDSHF : coefficient for activation of
$                              WAM4 dissipation for unsaturated (SDSLF) and
$                               saturated (SDSHF) parts of the spectrum
$                           SDSC2    : Saturation dissipation coefficient
$                           SDSC4    : Value of B0=B/Br for wich Sds is zero
$                           SDSBR    : Threshold Br for saturation
$                           SDSP     : power of (B/Br-B0) in Sds
$                           SDSBR2   : Threshold Br2 for the separation of
$                             WAM4 dissipation in saturated and non-saturated
$                           SDSC5 : coefficient for turbulence dissipation
$                           SDSC6 : Weight for the istropic part of Sds_SAT
$                           SDSDTH: Angular half-width for integration of B
$
$
$ Bottom friction  - - - - - - - - - - - - - - - - - - - - - - - - - -
$   JONSWAP             : Namelist SBT1
$                           GAMMA   : As it says.
$
$
$ Surf breaking  - - - - - - - - - - - - - - - - - - - - - - - - - - -
$   Battjes and Janssen : Namelist SDB1
$                           BJALFA  :
$                           BJGAM   :
$                           BJFLAG  :
$
$ Propagation schemes ------------------------------------------------ $
$   First order         : Namelist PRO1
$                           CFLTM  : Maximum CFL number for refraction.
$
$   UQ with diffusion   : Namelist PRO2
$                           CFLTM  : Maximum CFL number for refraction.
$                           DTIME  : Swell age (s) in garden sprinkler
$                                    correction. If 0., all diffusion
$                                    switched off. If small non-zero
$                                    (DEFAULT !!!) only wave growth
$                                    diffusion.
$                           LATMIN : Maximum latitude used in calc. of
$                                    strength of diffusion for prop.
$
$   UQ with averaging   : Namelist PRO3
$                           CFLTM  : Maximum CFL number for refraction.
$                           WDTHCG : Tuning factor propag. direction.
$                           WDTHTH : Tuning factor normal direction.
$
$ Miscellaneous ------------------------------------------------------ $
$   Misc. parameters    : Namelist MISC
$                           CICE0  : Ice concentration cut-off.
$                           CICEN  : Ice concentration cut-off.
$                           PMOVE  : Power p in GSE aleviation for
$                                    moving grids in Eq. (D.4).
$                           XSEED  : Xseed in seeding alg. (!/SEED).
$                           FLAGTR : Indicating presence and type of
$                                    subgrid information :
$                                     0 : No subgrid information.
$                                     1 : Transparancies at cell boun-
$                                         daries between grid points.
$                                     2 : Transp. at cell centers.
$                                     3 : Like 1 with cont. ice.
$                                     4 : Like 2 with cont. ice.
$                           XP, XR, XFILT
$                                    Xp, Xr and Xf for the dynamic
$                                    integration scheme.
$                           IHMAX  : Number of discrete levels in part.
$                           HSPMIN : Minimum Hs in partitioning.
$                           WSM    : Wind speed multiplier in part.
$                           WSC    : Cut of wind sea fraction for
$                                    identifying wind sea in part.
$                           FLC    : Flag for combining wind seas in
$                                    partitioning.
$                           FMICHE : Constant in Miche limiter.
$
$ In the 'Out of the box' test setup we run with sub-grid obstacles
$ and with continuous ice treatment.

$ Define grid
$ Four records containing :
$  1 NX, NY. As the outer grid lines are always defined as land
$    points, the minimum size is 3x3
$  2 Grid increments SX, SY (degr.or m) and scaling (division) factor
$    If NX*SX = 360., latitudinal closure is applied
$  3 Coordinates of (1,1) (degr.) and scaling (division) factor
$  4 Limiting bottom depth (m) to discriminate between land and sea
$    points, minimum water depth (m) as allowed in model, unit number
$    of file with bottom depths, scale factor for bottom depths (mult.),
$    IDLA, IDFM, format for formatted read, FROM and filename.
$      IDLA : Layout indicator :
$                  1   : Read line-by-line bottom to top
$                  2   : Like 1, single read statement.
$                  3   : Read line-by-line top to bottom
$                  4   : Like 3, single read statement
$      IDFM : format indicator :
$                  1   : Free format
$                  2   : Fixed format with above format descriptor
$                  3   : Unformatted
$      FROM : file type parameter
$               'UNIT' : open file by unit number only.
$               'NAME' : open file by name and assign to unit.

Dans mon cas, j'étais en IDLA : 1 et IDFM : 2  '(20f8.2)'.
Je pense donc qu'il commencait par lire le coin inferieur gauche
(qui correspond au sud ouest de la zone).

La routine fortran : (MECO, NECO etant Nx, Ny)

    program conversion_topo_symphonie_ww3
    PARAMETER(MECO=160,NECO=118)
        DIMENSION H_Z(MECO,NECO),MORPHO_Z(MECO,NECO)
        OPEN(UNIT=3,FILE='bathycote_in.dat',RECL=10000)
        DO I=1,MECO
        READ(3,'(1000I1)')(MORPHO_Z(I,J),J=1,NECO)
        ENDDO
        DO I=1,MECO
        READ(3,*)(H_Z(I,J),J=1,NECO)
        ENDDO
! format ww3
        OPEN(UNIT=20,FILE='bathy_finisterre_ww3.ascii')
        do j=1,neco
        write(20,'(20f8.2)')(h_z(i,j),i=1,meco)
        enddo
        close(20)

        OPEN(UNIT=20,FILE='masque_finisterre_ww3.ascii')
        do j=1,neco
        write(20,'(1000I1)')(morpho_z(i,j),i=1,meco)
        enddo
        close(20)
        end
-----------------------------------------------------------------------------*/

using namespace std;

class ww3_namelist_t{
private :
public :
/* *-------------------------------------------
$ Stresses
$   TC 1996 with cap    : Namelist FLX3
$                           CDMAX  : Maximum allowed CD (cap)
$                           CTYPE  : Cap type :
$                                     0: Discontinuous (default).
$                                     1: Hyperbolic tangent.
*/
#define WW3_FLX3                    0

#define WW3_FLX3_CDLMAX_LABEL       "CDMAX"
#define WW3_FLX3_CTYPE_LABEL        "CTYPE"

#define WW3_FLX3_CAP_DISCONTINUOUS  0
#define WW3_FLX3_CAP_HYPERBOLIC_TGT 1

  std::map<std::string, int> FLX3;

/* *-------------------------------------------
$ Linear input
$   Cavaleri and M-R    : Namelist SLN1
$                           CLIN   : Proportionality constant.
$                           RFPM   : Factor for fPM in filter.
$                           RFHF   : Factor for fh in filter.
*/
#define WW3_SLN1                    1

#define WW3_SLN1_CLIN_LABEL         "CLIN"
#define WW3_SLN1_RFPM_LABEL         "RFPM"
#define WW3_SLN1_RFHF_LABEL         "RFHF"

  std::map<std::string, int> SLN1;

/** to be continued ... */

/* *---------------------------------------------
  Exponential input
*/

/* *-------------------------------------------
$   WAM-3               : Namelist SIN1
$                           CINP   : Proportionality constant.
*/

/** to be continued ... */

/* *-------------------------------------------
$   Tolman and Chalikov : Namelist SIN2
$                           ZWND   : Height of wind (m).
$                           SWELLF : swell factor in (n.nn).
$                           STABSH, STABOF, CNEG, CPOS, FNEG :
$                                    c0, ST0, c1, c2 and f1 in . (n.nn)
$                                    through (2.65) for definition of
$                                    effective wind speed (!/STAB2).
*/

/** to be continued ... */

/* *-------------------------------------------
$   WAM4 and variants  : Namelist SIN3
$                           ZWND    : Height of wind (m).
$                           ALPHA0  : minimum value of Charnock coefficient
$                           Z0MAX   : maximum value of air-side roughness z0
$                           BETAMAX : maximum value of wind-wave coupling
$                           SINTHP  : power of cosine in wind input
$                           ZALP    : wave age shift to account for gustiness
$                       TAUWSHELTER : sheltering of short waves to reduce u_star
$                         SWELLFPAR : choice of swell attenuation formulation
$                                    (1: TC 1996, 3: ACC 2008)
$                            SWELLF : swell attenuation factor
$     Extra parameters for SWELLFPAR=3 only
$                  SWELLF2, SWELLF3 : swell attenuation factors
$                           SWELLF4 : Threshold Reynolds number for ACC2008
$                  	    SWELLF5 : Relative viscous decay below threshold
$                             Z0RAT : roughness for oscil. flow / mean flow
$
*/
#define WW3_SIN3                    4
#define WW3_SIN3_LABEL              "&SIN3"

#define WW3_SIN3_ZWND_LABEL         "ZWND"
#define WW3_SIN3_ALPHA0_LABEL       "ALPHA0"

  std::map<std::string, int> SIN3;

/** to be continued ... */

  ww3_namelist_t() {
    }

  void destroy() {
    }
};

class ww3_io_t{
private :
public :
/* *-------------------------------------------
  gridname */
  char gridname[30];
/* *-------------------------------------------
  frequencies and directions */
  double f0,df,nf,nd,doffset;
/* *-------------------------------------------
  time steps */
  double dt,cfl_xy,cfl_kt,dtsource_min;
/* *-------------------------------------------
  model flags */
  bool FLDRY, FLCX, FLCY ,FLCTH, FLCK, FLSOU;
/* *-------------------------------------------
  Stresses */
  double CDMAX;
  int    CTYPE;
/* *-------------------------------------------
  Linear input */
  double CLIN,RFPM,RFHF;

/** to be continued ... */

/* *-------------------------------------------
  Exponential input  */

// WAM-3 : Namelist SIN1

/** to be continued ... */

// Tolman and Chalikov : Namelist SIN2

/** to be continued ... */

// WAM4 and variants  : Namelist SIN3
  double ZWND,ALPHA0,Z0MAX,BETAMAX,SINTHP,ZALP,TAUWSHELTER;
  int    SWELLFPAR;
  double SWELLF;

/** to be continued ... */


/* *-------------------------------------------
  model grid*/
  size_t nx,ny;
  double SX,SY;
  double LON0,LAT0,scale;
  int    IDLA,IDFM;
  char   FROM[4];

  ww3_io_t() {
    nx=ny =0;
    }

  void destroy() {
    }
};

#define WW3_GRIDNAME_LINE     0
#define WW3_FREQUENCIES_LINE  1
#define WW3_TIME_STEPS_LINE   2
#define WW3_MODEL_FLAGS_LINE  3

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

//void ww3_parse_notebook_grid (const char *filename, ww3_io_t header)
void ww3_parse_notebook_grid (const char *filename)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  std::map<std::string, string> vars;
  std::string line;
  std::string t1, t2;
  std::ifstream file;
  std::vector<std::string> s1, s2;
  std::stringstream st1;
  unsigned int line_counter;
  ww3_io_t header;
/*-----------------------------------------------------------------------------
  open file */
  file.open(filename);

  line_counter=0;

  while(std::getline(file, line)) {
/*-----------------------------------------------------------------------------
    skip comment lines */
    if(line[0]=='$') continue;
    switch (line_counter) {
      case WW3_GRIDNAME_LINE:
//        header.gridname="\0";
        line_counter++;
        break;
      case WW3_FREQUENCIES_LINE:
        line_counter++;
        break;
      case WW3_TIME_STEPS_LINE:
        line_counter++;
        break;
      case WW3_MODEL_FLAGS_LINE:
        line_counter++;
        break;
      default:
        size_t pos;
        pos=line.find(WW3_SIN3_LABEL);
        if(pos!=-1) {
          }
        line_counter++;
        break;
      }
    std::string vname;
    int vvalue;
/* *----------------------------------------------------------------------------
    parse the line */
//    vars[vname] = vvalue;
    std::cout << line << std::endl;
    }
  file.close();
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void ww3_ascii_printheader (ww3_io_t header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  printf("#ascii_printheader \n");
  printf(" nx      : %d \n", header.nx);
  printf(" ny      : %d \n", header.ny);

  printf(" xmin    : %f \n", header.LON0);
  printf(" ymin    : %f \n", header.LAT0);

  printf(" dx      : %f \n", header.SX);
  printf(" dy      : %f \n", header.SY);
}


extern void ascii_skipline(FILE *in);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ww3_ascii_readheader (FILE *in,ww3_io_t *header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems,status=0;
  char c;

  nitems=fscanf(in,"%d %d"      , &header->nx,&header->ny);
  nitems=fscanf(in,"%lf %lf"    , &header->SX,&header->SY);
  nitems=fscanf(in,"%lf %lf %lf", &header->LON0,&header->LAT0,&header->scale);
  nitems=fscanf(in,"%d %d %s"   , &header->IDLA,&header->IDFM,header->FROM);

  ascii_skipline(in);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

   int ww3_ascii_writeheader (FILE *in,ww3_io_t header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nitems,status=0;
  nitems=fprintf(in,"%d %d"      , header.nx,header.ny);
  nitems=fprintf(in,"%lf %lf"    , header.SX,header.SY);
  nitems=fprintf(in,"%lf %lf %lf", header.LON0,header.LAT0,header.scale);
  nitems=fprintf(in,"%d %d %s"   , header.IDLA,header.IDFM,header.FROM);
  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

   int ww3_header_to_grid(ww3_io_t header, grid_t *grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;

  grid->nx=header.nx;
  grid->ny=header.ny;

  grid->xmin=header.LON0;
  grid->ymin=header.LAT0;

  grid->dx=header.SX;
  grid->dy=header.SY;

  grid->xmax=grid->xmin+(grid->nx-1)*grid->dx;
  grid->ymax=grid->ymin+(grid->ny-1)*grid->dy;

  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

   int ww3_grid_to_header(grid_t grid, ww3_io_t *header)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;

  header->nx=grid.nx;
  header->ny=grid.ny;

  header->LON0=grid.xmin;
  header->LAT0=grid.ymin;

  header->SX=grid.dx;
  header->SY=grid.dy;

  status=0;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int www3_ascii_loadgrid(char* filename,grid_t *grid, int status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int option;
  size_t offset, nread, size;
  float time;
  ww3_io_t header;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    status=-1;
    return(status);
    }

  status=ww3_ascii_readheader(in,&header);
  status=ww3_header_to_grid(header,grid);
  if (status !=0) {
    status=-1;
    return(status);
    }

  fclose(in);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int www3_ascii_readr1_01 (FILE *in,grid_t grid, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,nitems,status=0;
  char c;

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      nitems=fscanf(in,"%f",&(buffer[m]));
      if(nitems!=1) {
        printf("anomaly!!!\n");
        }
      }
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int www3_ascii_readr1_02 (FILE *in,grid_t grid, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,nitems,status=0;
  char c;

  for (j=grid.ny-1;j>-1;j--) {
    for (i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      nitems=fscanf(in,"%f",&(buffer[m]));
      if(nitems!=1) {
        printf("anomaly!!!\n");
        }
      }
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int www3_ascii_loadr1(char* filename,grid_t *grid, float *buffer, int status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int option;
  size_t offset, nread, size;
  float time;
  ww3_io_t header;

  status=0;

  in=fopen(filename,"rb");
  if (in == NULL) {
    status=-1;
    return(status);
    }

  status=ww3_ascii_readheader(in,&header);
  status=ww3_header_to_grid(header,grid);

  switch(header.IDLA) {
    case 1:
      status=www3_ascii_readr1_01 (in, *grid, buffer);
      break;
    case 2:
      status=www3_ascii_readr1_02 (in, *grid, buffer);
      break;
    case 3:
      break;
    case 4:
      break;
    default:
      break;
    }

  if (status !=0) {
    status=-1;
    return(status);
    }

  fclose(in);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ww3_ascii_saver1_01 (FILE *in,grid_t grid, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,nitems,status=0;
  char c;

  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      nitems=fscanf(in,"%f",&(buffer[m]));
      if(nitems!=1) {
        printf("anomaly!!!\n");
        }
      }
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ww3_ascii_saver1_02 (FILE *in,grid_t grid, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,nitems,status=0;
  char c;

  for (j=grid.ny-1;j>-1;j--) {
    for (i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      nitems=fscanf(in,"%f",&(buffer[m]));
      if(nitems!=1) {
        printf("anomaly!!!\n");
        }
      }
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int www3_ascii_saver1(char* filename,grid_t grid, float *buffer, int status)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  FILE *in;
  int option;
  size_t offset, nread, size;
  float time;
  ww3_io_t header;

  status=0;

  in=fopen(filename,"w");
  if (in == NULL) {
    status=-1;
    return(status);
    }

  status=ww3_grid_to_header(grid,&header);

  header.IDLA=1;
  header.IDFM=1;
  sprintf(header.FROM,"%s","NAME");

  status=ww3_ascii_writeheader(in,header);

  switch(header.IDLA) {
    case 1:
      status=ww3_ascii_saver1_01 (in, grid, buffer);
      break;
    case 2:
      status=ww3_ascii_saver1_02 (in, grid, buffer);
      break;
    case 3:
      break;
    case 4:
      break;
    default:
      break;
    }

  if (status !=0) {
    status=-1;
    return(status);
    }

  fclose(in);
  return(status);
}
