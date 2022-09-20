
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

<!-- USE FIND AND REPLACE TO EDIT THIS LINE, SO THAT print_help IS ALSO UPDATED -->
\brief tidal loading

<!-- A LINK TO main() or print_help() WILL NOT LINK TO THE RIGHT SOURCE ! -->
See the main function for how this works
and the print_help function for how to use this.
*/
/*----------------------------------------------------------------------------*/

#define MAIN_SOURCE

#include "version-macros.def" //for VERSION and REVISION
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <unistd.h>

#include <stdint.h> //uint64_t
#include <sys/time.h> //gettimeofday

#include "tools-structures.h"

#include "functions.h" //safely includes omp.h
#include "maths.h"
#include "map.h"
#include "geo.h"
#include "tides.h"
#include "fe-proto.h"
#include "poc-netcdf-data.hpp"
#include "datastream.h"
#include "topo.h"
#include "constants.h"

#include "fe-classes.h"


struct timeval stv;///<start timeval, for progression
//These global variables are used to count the amount of time this programme took for various tasks
struct timeval before;
double ct=0.;// computation

string cmd;///<command line


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  grid_t get_zonegrid_shifted(char *zone)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{

  bool zone_initialised;
  grid_t grid=get_zonegrid(zone,&zone_initialised);
 
  if(zone_initialised==false){
    double resolution=atof(zone);
    if(resolution==0) TRAP_ERR_EXIT(-1,"%s is neither a zone nor a resolution\n",zone);
    map_set2Dgrid(&grid,resolution*.5,-90+resolution*.5,+360.0,+90.0,resolution);
    }
  
  return(grid);

}


/*----------------------------------------------------------------------------*/
//option variables
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

char *inputPath=NULL,*meshPath=NULL,*zone=NULL,*topoPath=NULL,
  *outputPath=NULL,*debugPath=NULL;
const char *elevAmplitudeName="Ha",*elevPhaseName="Hg";
char *uAmpName=0,*uPhaName=0,*vAmpName=0,*vPhaName=0;
complex<double> *uTide,*vTide;
double *depth;

bool lsaOnly=false;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void printTotalMassWithFPE(complex<double> totalMass,complex<double> massFPE){
  printf("%g +/- %g kg, %g deg",abs(totalMass),sqrt((real(massFPE)+imag(massFPE))*.5),-arg(totalMass)*r2d);
  }
void printTotalMassWithFPE(double totalMass,double massFPE){
  printf("%g +/- %g kg",totalMass,massFPE);
  }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename To,typename Ti> To *massFromAreaAndElevation(int ni,double *area,Ti *elevation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// calculates complex mass from area and elevation
/**
\param ni number of input polygones
\param *area area of input polygones
\param *elevation average complex amplitude of the tide at the input polygone
\returns the complex masses
*/
/*----------------------------------------------------------------------------*/
{
  int ii;//input index
  
#if 0 /* some of this code does not compile with all versions of gcc */
  if(true){
    double *__restrict__ masses,*__restrict__ aa,*__restrict__ reta;
    masses=(double*)__builtin_assume_aligned(new double[ni],__BIGGEST_ALIGNMENT__);
    aa=(double*)__builtin_assume_aligned(area,__BIGGEST_ALIGNMENT__);
    reta=(double*)__builtin_assume_aligned(new double[ni],__BIGGEST_ALIGNMENT__);
    
    for(ii=0;ii<ni;ii++){
      reta[ii]=real(elevation[ii]);
      }
    
    STDERR_BASE_LINE("%p=%p*%p are hopefully %d-bytes-aligned. %d elements\n",masses,aa,reta,__BIGGEST_ALIGNMENT__,ni);
    
    extern char difftimer_unit_prefix;
    difftimer_unit_prefix='m';
    
    for(int run=0;run<2;run++){
      //if(run>0) initialize_OPENMP(2);
      
      difftimer_t difftimer(FILE_LINE_FUNC,9,run);
      
      #pragma omp parallel for if(run>0)
      for(ii=0;ii<ni;ii++){
        if(aa[ii]<=-NC_FILL_DOUBLE) TRAP_ERR_EXIT(ENOEXEC,"error\n");
        if(aa[ii]>=NC_FILL_DOUBLE) TRAP_ERR_EXIT(ENOEXEC,"error\n");
        masses[ii]=aa[ii]*reta[ii];
        }
      
      difftimer();
      
      #pragma omp parallel for if(run>0)
      for(ii=0;ii<ni;ii++){
        if(aa[ii]<=-NC_FILL_DOUBLE) TRAP_ERR_EXIT(ENOEXEC,"error\n");
        if(aa[ii]>=NC_FILL_DOUBLE) TRAP_ERR_EXIT(ENOEXEC,"error\n");
        masses[ii]=aa[ii]*reta[ii];
        }
      
      difftimer();
      
      #pragma omp parallel for if(run>0)
      for(ii=0;ii<ni;ii++){
        if(aa[ii]>=NC_FILL_DOUBLE) TRAP_ERR_EXIT(ENOEXEC,"error\n");
        masses[ii]=aa[ii]*reta[ii];
        }
      
      difftimer();
      
      #pragma omp parallel for if(run>0)
      for(ii=0;ii<ni;ii++){
        masses[ii]=aa[ii]*reta[ii];
        }
      
      difftimer();
      
      #pragma omp parallel for if(run>0)
      for(ii=0;ii<ni;ii++){
        masses[ii]=aa[ii]+reta[ii];
        }
      
      difftimer();
      
#define sn 2048
      
/*--------------------------------------------------------------------------------
      2nd fastest */
      
      #pragma omp parallel for if(run>0) schedule(static,sn)
      for(ii=0;ii<ni;ii++){
        masses[ii]=aa[ii]*reta[ii];
        }
      
      difftimer();
      
/*--------------------------------------------------------------------------------
      fastest */
      
      int si,si1;
      
      si1=ni/sn+1;
      
      #pragma omp parallel if(run>0)
      {
      int ii0,ii1;
      
      #pragma omp for
      for(si=0;si<si1;si++){
        ii0=si*sn;
        ii1=min(ni,ii0+sn);
        
        //#pragma simd /* Single Instruction Mutiple Data => ignored */
        for(ii=ii0;ii<ii1;ii++)
          masses[ii]=aa[ii]*reta[ii];
        
        }
      
      }
    
      difftimer();
      
/*--------------------------------------------------------------------------------
      2nd slowest */
      
      si1=ni/sn+1;
      
      #pragma omp parallel if(run>0)
      {
      double
        input1[sn] __attribute__((aligned(__BIGGEST_ALIGNMENT__))),
        input2[sn] __attribute__((aligned(__BIGGEST_ALIGNMENT__))),
        output[sn] __attribute__((aligned(__BIGGEST_ALIGNMENT__)));
      
      int ii0,ii1;
      int i,i1;
      
      #pragma omp for
      for(si=0;si<si1;si++){
        ii0=si*sn;
        ii1=min(ni,ii0+sn);
        i1=ii1-ii0;
        
        for(i=0,ii=ii0;ii<ii1;ii++,i++){
          input1[i]=aa[ii];
          input2[i]=reta[ii];
          }
        
        for(i=0;i<sn;i++)
          output[i]=input1[i]*input2[i];
        
        for(i=0,ii=ii0;ii<ii1;ii++,i++){
          masses[ii]=output[i];
          }
        
        }
      
      }
    
     
      }
    
    TRAP_ERR_EXIT(ENOEXEC,"testing\n");
    }
#endif
  
  To *masses,totalMass=0.,massFPE=0.;/*< masses and floating point error */
  double absMass,maxMass=-1.,totalArea=0.;
  const double
    FPE=pow(2,-52),
    EarthArea=4*M_PI*MeanEarthRadius2;
  //TRAP_ERR_EXIT(ENOEXEC,"just checkin double's FPE:%g\n",FPE);
  /** \bug 2012-05-22 Damien Allain : using a constant volumic mass */
  const double rho_water=1025.;/*< volumic mass of sea water */
  
  masses=new To[ni];
  
  STDOUT_BASE_LINE("Masses of %d points.\n",ni);
  gettimeofday(&before);
  
  for(ii=0;ii<ni;ii++){
    totalArea+=area[ii];
#define test_elevation_equal_1 0
#if test_elevation_equal_1
    elevation[ii]=1.;
#endif
    ///using \$ w=a \eta \rho \$
    masses[ii]=area[ii]*elevation[ii]*rho_water;
    if(isnan(masses[ii])){
      STDOUT_BASE_LINE("");printc("g",masses[ii]);printf("=%g*",area[ii]);printc("g",elevation[ii]);printf("*%g\n",rho_water);
      TRAP_ERR_EXIT(ENOEXEC,"NAN found at %d\n",ii);
      }
    totalMass+=masses[ii];
    massFPE+=square(totalMass*FPE);
    absMass=abs(masses[ii]);
    updatemax(&maxMass,absMass);
    }
  
  const char *realInputPath=0;
  if(inputPath!=0){
    realInputPath=inputPath;
    }
  else if(zone!=0){
    realInputPath=zone;
    }
  
  STDOUT_BASE_LINE("Took %gs.\n",difftime(before));
  
  if(isnan(totalMass) or inputPath==0){
    ii=printf("Mass report for %s\n",realInputPath);
    for(;ii>1;ii--)printf("=");printf("\n");
    printf("All |masses|<= %g kg.\nTotal (that should ideally be 0kg): ",maxMass);
    printTotalMassWithFPE(totalMass,massFPE);
    printf(" over %g%% of Earth\n",100.*totalArea/EarthArea);
    
    if(isnan(totalMass))
      TRAP_ERR_EXIT(ENOEXEC,"See above.\n");
    }
  
  return masses;
}


/*----------------------------------------------------------------------------*/
///table of Green functions
/**
The columns are :
- angular distance (deg)
- radial deformation or displacement (m*1e12*6371e6*rad/kg)
- effet on vertical component of gravity
- effet on the gravitational potential
*/
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

#define GREEN_DEFORMATION 1
#define GREEN_DISPLACEMENT 2

//#define GREEN_NCOLS 3
//#define GREEN_NCOLS 4
#define GREEN_NCOLS 8

const double Green[][GREEN_NCOLS]={
#if GREEN_NCOLS == 3

#define POT_UP_TO_LSA 1
#define POT_OR_LSA_COLUMN 2
#define HORIZONTAL_COLUMN -1
#define GREEN_MODE GREEN_DISPLACEMENT

#include "/data1/LSA/FES2014/Boy/sal_CMEA.h"

#elif GREEN_NCOLS == 4

#define POT_UP_TO_LSA 1
#define POT_OR_LSA_COLUMN 3

#if 0
/**
Hard-coded from \c PAGIA
(Pagiatakis 1990 http://dx.doi.org/10.1111/j.1365-246X.1990.tb01790.x, AppB)
see Lubes PhD (1998, eqs. 11, 12 and 13 and Tabs. 3 and 4)
but automatically generated with :
\code
sed -re 's/( *[^ ]+)( +[^ ]+)( +[^ ]+)( +[^ ]+)/  {\1,\2,\3,\4}/;1 s/^/[][4]={\n/;$ {s/$/\n  };/;q};s/$/,/' PAGIA
\endcode
See also Longmann (1963 http:dx.doi.org/10.1111/j.1365-246X.1990.tb01790.x) or maybe Farrell (1972 http:dx.doi.org/10.1029/RG010i003p00761, TabA4 and TabA5)
*/
#define GREEN_MODE GREEN_DEFORMATION
#define HORIZONTAL_COLUMN -1
  {  0.0001, -42.6038,  -98.875,     -.00033},
  {  0.001,  -42.3426,  -98.361,     -.00248},
  {  0.01,   -39.7564,  -93.260,     -.01645},
  {  0.02,   -37.8453,  -87.757,     -.02789},
  {  0.03,   -35.8124,  -82.552,     -.03749},
  {  0.04,   -33.2968,  -77.732,     -.04492},
  {  0.06,   -29.2265,  -69.374,     -.06050},
  {  0.08,   -26.6460,  -62.674,     -.07307},
  {  0.10,   -24.3677,  -57.490,     -.08437},
  {  0.16,   -20.2121,  -49.808,     -.11408},
  {  0.20,   -19.2689,  -46.668,     -.13206},
  {  0.25,   -18.5960,  -45.043,     -.15338},
  {  0.30,   -17.5849,  -41.839,     -.17367},
  {  0.40,   -15.6612,  -34.798,     -.21140},
  {  0.50,   -14.8092,  -31.531,     -.24560},
  {  0.60,   -14.8010,  -31.690,     -.27664},
  {  0.80,   -14.6419,  -31.763,     -.33079},
  {  1.0,    -13.8089,  -29.756,     -.37654},
  {  1.2,    -12.9244,  -27.805,     -.41580},
  {  1.6,    -11.3105,  -24.374,     -.47997},
  {  2.0,    -10.1118,  -22.005,     -.53050},
  {  2.5,     -8.7091,  -19.033,     -.58081},
  {  3.0,     -7.7017,  -16.953,     -.62140},
  {  4.0,     -6.1931,  -13.671,     -.68449},
  {  5.0,     -5.2462,  -11.534,     -.73311},
  {  6.0,     -4.6370,  -10.081,     -.77306},
  {  7.0,     -4.2289,   -9.045,     -.80700},
  {  8.0,     -3.9542,   -8.309,     -.83615},
  {  9.0,     -3.7559,   -7.747,     -.86108},
  { 10.0,     -3.6083,   -7.308,     -.88205},
  { 12.0,     -3.3904,   -6.667,     -.91254},
  { 16.0,     -3.0224,   -5.714,     -.92975},
  { 20.0,     -2.6540,   -4.946,     -.89276},
  { 25.0,     -2.1288,   -4.018,     -.77938},
  { 30.0,     -1.5363,   -3.059,     -.60346},
  { 40.0,     -0.2903,   -1.468,     -.12340},
  { 50.0,      0.8640,   -0.456,      .419871},
  { 60.0,      1.6848,    0.159,      .90526},
  { 70.0,      2.0930,    0.512,     1.24070},
  { 80.0,      2.0506,    0.181,     1.37116},
  { 90.0,      1.6135,   -0.427,     1.27629},
  {100.0,      0.8960,   -0.890,      .96372},
  {110.0,     -0.0494,   -1.528,      .46169},
  {120.0,     -1.1269,   -2.215,     -.18716},
  {130.0,     -2.2614,   -2.768,     -.93077},
  {140.0,     -3.4025,   -3.272,    -1.71134},
  {150.0,     -4.4795,   -3.749,    -2.46841},
  {160.0,     -5.4136,   -4.133,    -3.14143},
  {170.0,     -6.1454,   -4.412,    -3.67235},
  {180.0,     -6.6420,   -4.660,    -4.00830}
#else
#define GREEN_MODE GREEN_DEFORMATION
#define HORIZONTAL_COLUMN -1
#include "REF_6371_loading_green_functions.h"
#endif

#elif GREEN_NCOLS == 8

#define POT_UP_TO_LSA -1
#define POT_OR_LSA_COLUMN 5

#define GREEN_MODE GREEN_DISPLACEMENT
#define HORIZONTAL_COLUMN 2 /* twice as slow when enabled */
#include "REF_6371_loading_green_functions-CoM.h"

#else
#error See lines above
#endif

#if HORIZONTAL_COLUMN <= 0
#warning NOT COMPUTING HORIZONTAL DEFORMATION !
#endif
  };
const int Greenn=sizeof(Green)/sizeof(*Green);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


/*----------------------------------------------------------------------------*/
//sqrt table
/*----------------------------------------------------------------------------*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
#define TABULATE_sqrt 0
#if TABULATE_sqrt
class tabledSqrt_t{
  double *table,idx;
  int n;
  
public:
#define TEST_tabledSqrt_t 0
  
  tabledSqrt_t(int n0=-1){
    if(n0<=0)n0=0x300001;
    n=n0;
    table=new double[n];
    idx=n/.75;
    const double dx=.75/n;
    for(int i=0;i<n;i++){
      double x=.25+i*dx;
      table[i]=sqrt(x);
      }

#if TEST_tabledSqrt_t
    const double kw=cbrt(2.);
    for(double x=.125;x<17.;x*=kw)
      operator()(x,1);
#endif
    }
  
  ~tabledSqrt_t(){
    delete[]table;
    }
  
  double operator()(double x
#if TEST_tabledSqrt_t
                    ,int verbose=0
#endif
                    ){
    /* NOTE: this is NOT faster */
    /* TODO: time budget to find out why... */
    int e;
    double fx=frexp(x,&e);
    if(e&1){fx*=.5;e++;}
    
    int i=(fx-.25)*idx;
    const double
      fr=table[i],
      r=ldexp(fr,e/2);
    
#if TEST_tabledSqrt_t
    if(verbose)STDOUT_BASE_LINE("tsqrt(%g)=%g\n",x,r);
#endif
    
    return r;
    }
  };

tabledSqrt_t tsqrt;
#else
#define tsqrt sqrt
#endif
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/


typedef union{
  double d;
  uint64_t i;
  } flog2_t;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int TEST_flog2_t()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  flog2_t t;
  const double
    fs[]={1,1.125,1.25,1.5};
  const int
    n=sizeof(fs)/sizeof(*fs);
  
  for(double s=-1;s<=1;s+=2.)
    for(int e=-2;e<=2;e++)
      for(int i=0;i<n;i++){
        //TRAP_ERR_EXIT(ENOEXEC,"s=%g;e=%g;i=%d\n",s,e,i);
        t.d=s*exp2(e)*fs[i];
        printf("%17g:%lx\n",t.d,t.i);
        }
  
  TRAP_ERR_EXIT(ENOEXEC,"%s is for testing\n",__func__);
  
  return 0;
}
//int test_flog2_t=TEST_flog2_t();


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

inline double flog2(double x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// fast approximation of log2, with an offset and an error < 0.09
/*----------------------------------------------------------------------------*/
{
#define CUSTOM_flog2 1
#if CUSTOM_flog2 == 1
/* MASK FORM
  1.6 times faster than frexp form
  with an offset of 1024 */
  flog2_t f,fx;
  f.d=x;
  
  static const int64_t
    frac_mask=0x000fffffffffffff,
    exp_mask =0x7ff0000000000000,
    //sign_mask=0x8000000000000000,/* unused */
    get_mask =0x3ff0000000000000;
  
  const int16_t
    e=(f.i&exp_mask)>>52;
  
  fx.i=(f.i&frac_mask)|get_mask;
  
  /* NOTE: ignoring sign bit, which, if set, would require a NAN to be returned */
  
  return e+fx.d;
#elif CUSTOM_flog2 == 0
/* FREXP FORM
  3 times faster than slow form
  with an offset of 2 */
  int e;
  const double
    fx=frexp(x,&e);
  
  return e+2*fx;
#elif CUSTOM_flog2 == -1
/* slow but exact form */
  return log2(x);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double fexp2(double x)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// reverse of flog2, not really optimised because this is not important
/*----------------------------------------------------------------------------*/
{
#if CUSTOM_flog2 >= 0
  double e,
    fx=modf(x,&e);
  
#if CUSTOM_flog2 == 1
  e-=1023;
#elif CUSTOM_flog2 == 0
  e--;
#endif
  
  fx++;
  fx/=2;
  
  return ldexp(fx,e);
#elif CUSTOM_flog2 == -1
  return exp2(x);
#endif
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int TEST_flog2()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const char *outPath=
#if CUSTOM_flog2 == 1
    "MASK_"
#elif CUSTOM_flog2 == 0
    "FREXP_"
#elif CUSTOM_flog2 == -1
    "SLOW_"
#endif
    "flog2.test";
  
  FILE
    *F=fopen(outPath,"w");
  if(F==NULL) TRAP_ERR_EXIT(errno,"failed to open %s (%d %m)\n",outPath,errno);
  
  STDOUT_BASE_LINE_FUNC("Ouput is %s\n",outPath);
  
  const double
    fs[]={1,1.125,1.25,1.5};
  const int
    n=sizeof(fs)/sizeof(*fs);
  
  for(int e=-2;e<=2;e++)
    for(int i=0;i<n;i++){
      //TRAP_ERR_EXIT(ENOEXEC,"s=%g;e=%g;i=%d\n",s,e,i);
      double
        x0=exp2(e)*fs[i],
        lx=flog2(x0),
        x1=fexp2(lx);
      fprintf(F,"%17g %17g %17g\n",x0,lx,x1);
      }
  
  fclose(F);
  
  TRAP_ERR_EXIT(ENOEXEC,"%s is for testing\n",__func__);
  
  return 0;
}
//int test=TEST_flog2();


#define PROPAGATION_DELAY 0


typedef
#if PROPAGATION_DELAY
  complex<double>
#else
  double
#endif
  radial_deformation_t;


typedef struct{
  radial_deformation_t
#if HORIZONTAL_COLUMN > 0
    h,/*< horizontal deformations */
#endif
    r;/*< radial deformations */
  double gp;/*< gravitational potential */
  } GreenTableRow_t;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T getCADTable(T *CumSum,const T & dCumSum,double d2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *CumSum+=dCumSum;
  return (2./d2)* *CumSum;
}


#if HORIZONTAL_COLUMN > 0
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> inline void h2en(const vector3_t & euv, const vector3_t & nuv, const vector3_t & imuv, const T &hord, T *east, T *north)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  convert horizontal displacement to eastern and northern components

  euv : eastern unitary vector
  nuv : northern unitary vector
  imuv : unitary vector of the input mass
  hord : horizontal displacement
  *eastio : eastern component. Incremented.
  *nortio : nortern component. Incremented.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double e,n,k;
  
  e=imuv*euv;
  n=imuv*nuv;
  k=1/hypot(e,n);
  e*=k;
  n*=k;
  
  *east+=hord*e;
  *north+=hord*n;
}
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> inline void lsa_increment(const GreenTableRow_t *table, int it, const T & mass, const vector3_t & euv, const vector3_t & nuv, const vector3_t & imuv, T *up, T *east, T *north, T *pot, T *lsa)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const GreenTableRow_t *tablei=&table[it];/* 13.3s */
  
  if(not lsaOnly){
    *up+=mass*tablei->r;/* 42s */
    
#if HORIZONTAL_COLUMN > 0
    T hord;
    hord=mass*tablei->h;/* 25s */
    h2en(euv,nuv,imuv,hord,east,north);/* 137s */
#endif
    
#if POT_UP_TO_LSA == 1
    *pot +=mass*tablei->gp;/* 32s */
#endif
    }
  
#if POT_UP_TO_LSA == -1
  *lsa +=mass*tablei->gp;
#endif
}


#if PROPAGATION_DELAY
//#include "tools-structures.h"
#include "tides.def"
#endif


#define CLUSTERING_ENABLED 1

/* Define to 1 to test a (rubish) fast gradient method */
#define FAST_GRAD 0


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class V,class M,typename T> double loading(int ni,const double *lai,const double *loi,const T *mass,int no,const double *lao,const double *loo,T *deformation,T *displacement,T *east,T *north,T *pot,T *lsa,mesh_t *mesh=NULL,grid_t *grid=NULL)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Computes radial deformation of the Earth by a tidal wave and gravitational potential of the tidal wave
/**
\param ni number of input points
\param *lai latitude in degrees of input points
\param *loi longitude in degrees of input points
\param *mass complex mass of the input points
\param no number of output points
\param *lao latitude in degrees of output points
\param *loo longitude in degrees of output points
\param *deformation complex amplitude of the radial deformation at the output points. Can be NULL.
\param *displacement complex amplitude of the radial displacement at the output points
\param *east complex amplitude of the eastern displacement at the output points: see note below
\param *north complex amplitude of the northern displacement at the output points: see note below
\param *pot complex amplitude of the gravitational potential at the output points
\param *lsa lsa[]=pot[]-displacement[]
\param *mesh mesh, if relevant.
\param *grid grid, if relevant.
\param shortestEgdeL2 in rad^2, if the mesh is not given
\returns core computation time

\note \c north and \c east are only computed if the horizontal Green function is activated
*/
/*----------------------------------------------------------------------------*/
{
  int ii,io;//input and output indexes
  static vector3_t *ui=0;//UNITARY coordinates of inputs
  V momentumDisplacement;
  #warning cvector3_t will have to be replaced by V in line below
  cvector3_t angularMomentum;
  M inertiaTensor;
  static vector3_t *uo=0;//UNITARY coordinates of outputs
  double cLat=NAN,sLat=NAN,cLon=NAN,sLon=NAN;//sines and cosines of latitude and longitude
  
/*---------------------------------------------------------------------*//**<h1>
  computes unitary x, y and z </h1>*/
  gettimeofday(&before);
  bool doUi=false;
  if(ui==0){
    STDOUT_BASE_LINE("Coordinates:\n");
    STDOUT_BASE_LINE("- unitary xi, yi and zi (and momenta)...\n");
    ui=new vector3_t[ni];
    
    doUi=true;
    }
  
  for(ii=0;ii<ni;ii++){
    
    if(doUi or uTide!=0 and vTide!=0){
      cLat=cos(lai[ii]*d2r);
      sLat=sin(lai[ii]*d2r);
      cLon=cos(loi[ii]*d2r);
      sLon=sin(loi[ii]*d2r);
      }
    
    vector3_t *uii=&ui[ii];
    
    if(doUi)
      uii->init(cLat*cLon,cLat*sLon,sLat);
    
    momentumDisplacement+=(*uii)*mass[ii];
    
    if(uTide!=0 and vTide!=0){
      const vector3_t
        uEast (-sLon     , cLon     ,0.),
        uNorth(-sLat*cLon,-sLat*sLon,cLat);
      /* What follows could be optimised because the ENZ base is direct orthonormal,
        making unneeded the vectorial product in the optimised version. */
      const cvector3_t
        waterColumMomentum=uTide[ii]*uEast+vTide[ii]*uNorth;
      angularMomentum+=*uii & waterColumMomentum;
      }
    
    /* See https://fr.wikipedia.org/wiki/Moment_d%27inertie#Mise_en_%C3%A9vidence_et_d%C3%A9finition_du_tenseur_d%27inertie */
    
    inertiaTensor.c[0][0]+=mass[ii]*(square(uii->y)+square(uii->z));
    inertiaTensor.c[0][1]+=mass[ii]*uii->x*uii->y;
    inertiaTensor.c[0][2]+=mass[ii]*uii->x*uii->z;
    
    inertiaTensor.c[1][1]+=mass[ii]*(square(uii->x)+square(uii->z));
    inertiaTensor.c[1][2]+=mass[ii]*uii->y*uii->z;
    
    inertiaTensor.c[2][2]+=mass[ii]*(square(uii->x)+square(uii->y));
    }
  
  inertiaTensor.c[1][0]=inertiaTensor.c[0][1];
  inertiaTensor.c[2][0]=inertiaTensor.c[0][2];
  inertiaTensor.c[2][1]=inertiaTensor.c[1][2];
  
  const double EarthMass=5.9722e+24;/* from http://dx.doi.org/10.1007/s10569-011-9352-4 Tab1 */
  momentumDisplacement*=MeanEarthRadius/EarthMass;
  /* momentumDisplacement is now
  the position of the center of mass with reference to the center of figure */
  
  angularMomentum*=MeanEarthRadius;
  inertiaTensor  *=square(MeanEarthRadius);
  
  if(outputPath==0){
    STDOUT_BASE_LINE("\n"
      "Momentum displacement\n"
      "=====================\n"
      "Displacement vector components in m are:\n");
    printf("  * x: ");printc("g",momentumDisplacement.x);printf("\n");
    printf("  * y: ");printc("g",momentumDisplacement.y);printf("\n");
    printf("  * z: ");printc("g",momentumDisplacement.z);printf("\n");
    
    double major,a;
    ellipse_Madmp3D(momentumDisplacement.x,momentumDisplacement.y,momentumDisplacement.z,&major,&a);
    printf("so the maximum displacement is of %g m %g deg.\n\n",major,a);
    
    if(uTide!=0 and vTide!=0){
      STDOUT_BASE_LINE("\n"
        "Relative angular momentum\n"
        "=========================\n"
        "In kg.m2.s-1:\n");
      printf("  * x: ");printc("g",angularMomentum.x);printf("\n");
      printf("  * y: ");printc("g",angularMomentum.y);printf("\n");
      printf("  * z: ");printc("g",angularMomentum.z);printf("\n");
      }
    
    STDOUT_BASE_LINE("\n"
      "Inertia tensor\n"
      "==============\n"
      "In kg.m2:\n");
    
    int i,j;/*< tensor indexes */
    int col;/*< column */
    for(j=0;j<3;j++){
      col=0;
      for(i=0;i<3;i++){
        col+=printf("  ");
        col+=printc("12e",inertiaTensor.c[j][i]);
        }
      printf("\n");
      }
    
    TRAP_ERR_EXIT(0,"No output path given: just checking masses and momenta.\n");
    }
  
  if(uo==0){
    STDOUT_BASE_LINE("- unitary xo, yo and zo...\n");
    uo=new vector3_t[no];
    for(io=0;io<no;io++){
      cLat=cos(lao[io]*d2r);
      vector3_t *uoi=&uo[io];
      uoi->x=cLat*cos(loo[io]*d2r);
      uoi->y=cLat*sin(loo[io]*d2r);
      uoi->z=sin(lao[io]*d2r);
      }
    }
  
  for(io=0;io<no;io++){
    if(not lsaOnly){
      displacement[io]=NC_FILL_DOUBLE;
      if(deformation!=0)
        deformation[io]=NC_FILL_DOUBLE;
#if HORIZONTAL_COLUMN > 0
      east[io]=NC_FILL_DOUBLE;
      north[io]=NC_FILL_DOUBLE;
#endif
      pot[io]=NC_FILL_DOUBLE;
      }
    lsa[io]=NC_FILL_DOUBLE;
    }
  STDOUT_BASE_LINE("Took %gs.\n",difftime(before));
  
/*----------------------------------------------------------------------------*/
  static double shortestEgdeL2=NAN;
  
  const double
    maxFlog2r2=flog2(4.);
  static double
    maxFracD2=4.;
  
  if(isnan(shortestEgdeL2)){
    if(mesh!=NULL){
      if(mesh->vertices[0].elmts==NULL){
        /* build element lists if they are missing */
        STDOUT_BASE_LINE("Element tables of vertices.\n");
        fe_vertex_element_tables(mesh);
        }
      
      /* find longest and shortest edges */
      range_t<double> l2R;/*< (rad^2) */
      vector2D_t shortestLoLa,longestLoLa;
      bool updateShortest,updateLongest;
      for(int ie=0;ie<mesh->nedges;ie++){
        edge_t *edge=&mesh->edges[ie];
        const int
          v0=edge->extremity[0],
          v1=edge->extremity[1];
        const double
          /* square of distance */
          //l2=edge->L;//FAILS
          l2=square(uo[v0]-uo[v1]);
        updateShortest=updatemin(&l2R.min,l2);
        updateLongest =updatemax(&l2R.max,l2);
        if(updateShortest or updateLongest){
          edge->lon=mesh->vertices[v0].lon;
          edge->lon+=degree_recale(mesh->vertices[v1].lon,edge->lon);
          edge->lon*=.5;
          edge->lat=(mesh->vertices[v0].lat+mesh->vertices[v1].lat)*.5;
          }
        if(updateShortest)
          shortestLoLa.init(edge->lon,edge->lat);
        if(updateLongest )
          longestLoLa .init(edge->lon,edge->lat);
        }
      
      /* NOTE: this will be over what flat-ish triangles will give */
      shortestEgdeL2=l2R.min/12;/* (cos(M_PI/3)/3)^2 = 1/12 */
      STDOUT_BASE_LINE("Shortest edge is at lo=%-.3f;la=%-.3f;zo=1e4 and is %grad=%gm so flog2(r*r) in [%g;%g](%g)\n",
        shortestLoLa.x,shortestLoLa.y,
        sqrt(l2R.min),sqrt(l2R.min)*MeanEarthRadius,
        flog2(l2R.min),maxFlog2r2,maxFlog2r2-flog2(l2R.min));
      
      maxFracD2=square(10)*l2R.max;
      STDOUT_BASE_LINE(" Longest edge is at lo=%-.3f;la=%-.3f;zo=1e2 and is %grad=%gm. *10=%grad=%gm\n",
        longestLoLa.x,longestLoLa.y,
        sqrt(l2R.max),sqrt(l2R.max)*MeanEarthRadius,
        sqrt(maxFracD2),sqrt(maxFracD2)*MeanEarthRadius);
      }
    if(grid!=0){
      shortestEgdeL2=square(grid->dy*d2r*grid->dx/360.);
      maxFracD2     =square(10*grid->dy*d2r);
      }
    STDERR_BASE_LINE(">=%g rad = %g deg : %d %d; frac if > %g rad = %g deg\n",
      sqrt(shortestEgdeL2),sqrt(shortestEgdeL2)*r2d,shortestEgdeL2>0,not(shortestEgdeL2>0),sqrt(maxFracD2),sqrt(maxFracD2)*r2d);
    }
  if(not(shortestEgdeL2>0)){
    TRAP_ERR_EXIT(ENOEXEC,"programming error\n");
    }
  
/*------------------------------------------------------------------------------
  tabulate Green function for radial deformation and for gravitational potential factor */
  int it;//table index
  static GreenTableRow_t *table=0;     /*< table of deformations + gravitational potential */
  static GreenTableRow_t *CADTable=0;  /*< table of cumulated averages multiplied by distance */
  const long
    /* We want 1% error, but no rounding error, so that is
    2^7=128 numbers per unit of flog2(r)~flog2(r**r)/2, or
    64 numbers per unit of flog2(r**r) . */
    tableInterpolation=64,
    tableLength=(maxFlog2r2-flog2(shortestEgdeL2))*tableInterpolation+1;
  
  if(tableLength<=0) TRAP_ERR_EXIT(ENOEXEC,"programming error: %d=(%g-flog2(%g)=%g)*%d+1\n",tableLength,maxFlog2r2,shortestEgdeL2,flog2(shortestEgdeL2),tableInterpolation);
  
  if(table==0 or CADTable==0){
    if(table!=0 or CADTable!=0) TRAP_ERR_EXIT(ENOEXEC,"programming error: %p!=0 or %p!=0\n",table,CADTable);
    
    STDOUT_BASE_LINE("Tabulating %d Green function values for load from %d values...\n",
      tableLength,Greenn);
    table   =new GreenTableRow_t[tableLength];
    CADTable=new GreenTableRow_t[tableLength];
    gettimeofday(&before);
    double angle,d[Greenn];//angle, Euclidian distance on the unitary sphere
    double r,rd[Greenn];//radial deformation
#if HORIZONTAL_COLUMN > 0
    double h,hd[Greenn];//horizontal deformation
#endif
    double gp,gpd[Greenn];//gravitational potential
    
    /// <h2>For this, it first converts the values given by ::Green</h2>
    
#if PROPAGATION_DELAY
    tidal_wave wave=wK1;
    wave.init();
    STDERR_BASE_LINE_FUNC("PROPAGATION_DELAY OF %s (%g deg/h)\n",wave.name,wave.omega);
    if(not isnormal(wave.omega))wexit(ENOEXEC);
    const double
      v=5e3,/* wave speed in solid Earth (m/s) */
      k=wave.omega*dph2rps/v,/* wave number in solid Earth (rad/m) */
      kMER=k*MeanEarthRadius;
#endif
    
    for(it=0;it<Greenn;it++){
      angle=Green[it][0]*d2r;//rad
      /* Euclidian distance on the unitary sphere \f$ d=2\sin\frac{\alpha}{2} \f$ */
      d[it]=2.*sin(angle*.5);
      /* converting the deformations to m.kg-1 */
#if GREEN_NCOLS == 3
      r=Green[it][1]/(-sin(angle)*MeanEarthRadius2/9.79764322);
#else
      r=Green[it][1]/(angle*MeanEarthRadius*1e12);
#endif
      /* and multiplying both */
      rd[it]=r*d[it];
#if HORIZONTAL_COLUMN > 0
      h=Green[it][HORIZONTAL_COLUMN]/(angle*MeanEarthRadius*1e12);
      hd[it]=h*d[it];
#endif
      /**- calculating the gravitational potential factor with \code /**/ // COMPILED CODE BELOW !!!
#if GREEN_NCOLS == 4
      /* PAGIA-style :
      g=9.81=G*M/R/R <=> 1/R*G/g=1/R*G/(G*M/R/R)=R/M */
      gpd[it]=(1+Green[it][3]*.149)*MeanEarthRadius/5.9736e24;
#else
      gp=Green[it][POT_OR_LSA_COLUMN]/(angle*MeanEarthRadius*1e12);
      gpd[it]=gp*d[it];
#endif
      /** \endcode */
      //STDERR_BASE_LINE("[%d] %g*%g=%g\n",it,r,d[it],rd[it]);//wexit(ENOEXEC);
      }
    /// <h2>It then interpolates.</h2>
    double flog2d2,d2,dit;
    int acc=-1,status;
    for(it=0;it<tableLength;it++){
      flog2d2=maxFlog2r2-(it+.5)/tableInterpolation;
      d2=fexp2(flog2d2);
      dit=sqrt(d2);
      
      status=map_interpolate1D(rd ,d,NAN,Greenn,dit,&table[it].r ,0,&acc);
      table[it].r/=dit;
#if PROPAGATION_DELAY
      table[it].r*=polar(1.,kMER*hypot(dit,0.1));
#endif
      if(status!=0 or isnan(table[it].r )) TRAP_ERR_EXIT(status,"%d %g %g\n",it,dit,d[0]);
      
#if HORIZONTAL_COLUMN > 0
      status=map_interpolate1D(hd ,d,NAN,Greenn,dit,&table[it].h ,0,&acc);
      table[it].h/=dit;
#if PROPAGATION_DELAY
    #error not coded yet
#endif
      if(status!=0 or isnan(table[it].h )) TRAP_ERR_EXIT(status,"%d %g %g\n",it,dit,d[0]);
#endif
      
      status=map_interpolate1D(gpd,d,NAN,Greenn,dit,&table[it].gp,0,&acc);
      table[it].gp/=dit;
      if(status!=0 or isnan(table[it].gp)) TRAP_ERR_EXIT(status,"%d %g %g\n",it,dit,d[0]);
      }
    
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

compute the double of the cumulated average Green function for radial deformation
because otherwise differences in the sizes of the triangles induce noise in the output.

The maths is as follows. A triangle is approximated to the portion of a pie, the deformation being calculated at the center of the arc.
Taking a portion of a pie of mass \f$ m \f$, of angle \f$ \alpha \f$ and of radius \f$ r \f$ gives its area \f$ A \f$:
\f[ A = \int_0^r{\alpha x \partial x} = 2^{-1}\alpha r^2 \f]
Taking the mass distribution \f$ p = m a^{-1} \f$:
\f[ m = A p = 2^{-1}\alpha r^2 p
\Leftrightarrow \alpha p = 2 m r^{-2}\f] and \f$ g(x)=G(x)x^{-1} \f$ the green functions gives the deformation \f$ D \f$:
\f[ D = \int_0^r{\alpha x p g(x) \partial x}
= \alpha p \int_0^r{G(x) \partial x}
= 2 m r^{-2} \int_0^r{G(x) \partial x}
= m r^{-1} \left[ 2 \frac{1}{r} \int_0^r{G(x) \partial x} \right] \f]
with the term in square brackets being the double of the cumulated average Green function for radial deformation.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
    FILE *f=NULL;
    //f=fopen("CADTable.dat","w");
    
    double dit0=0.,dd,
      GPDCumSum=0.;
    radial_deformation_t
#if HORIZONTAL_COLUMN > 0
      HDDCumSum=0.,
#endif
      RDDCumSum=0.;
    
    for(it=tableLength-1;it>=0;it--){
      flog2d2=maxFlog2r2-(it+.5)/tableInterpolation;
      d2=fexp2(flog2d2);
      dit=sqrt(d2);
      dd=dit-dit0;
      dit0=dit;
      
      CADTable[it].r =getCADTable(&RDDCumSum, table[it].r *dit*dd, d2);
#if HORIZONTAL_COLUMN > 0
      CADTable[it].h =getCADTable(&HDDCumSum, table[it].h *dit*dd, d2);
#endif
      CADTable[it].gp=getCADTable(&GPDCumSum, table[it].gp*dit*dd, d2);
      
      if(f!=NULL){
        fprintf(f,"%g %g %g %g %g",dit,
          table[it].r *dit,CADTable[it].r ,
          table[it].gp*dit,CADTable[it].gp);
#if HORIZONTAL_COLUMN > 0
        fprintf(f," %g %g",
          table[it].h *dit,CADTable[it].h );
#endif
        fprintf(f,"\n");
        }
      }
    
    if(f!=NULL) TRAP_ERR_EXIT(ENOEXEC,"testing\n");
    
    STDOUT_BASE_LINE("Took %g s.\n",difftime(before));
    }
  
/*------------------------------------------------------------------------------
  clustering masses */
  
#if CLUSTERING_ENABLED
  gettimeofday(&before);
  
// #warning comment this line and put value below back to 10.
  const double sf=10.;/*< security factor */
  /*for ca=1., we have, for the following values of sf:
      + 20: 5µm 345s
      + 10:13µm 242s
      +  5:80µm 202S
    */
  double ca;/*< clustering angle, in degrees.
    Should be of the order of (in radians) :
      + (4*pi/(sf*ni^.5))^.5 for reduced Gaussian grid
      + pi*(2/sf/(2*ni/pi)^.5)^.5 for regular grid
    but as a test on its closeness is made with UGs,
    this is about half this:
      + 0.5°: 400 s
      + 1.0°: 226.345 s
      + 1.5°: 226.219 s */
  ca=M_PI*sqrt(2/sf/sqrt(2*ni/M_PI))*r2d;
  int log2ca_x2=round(log2(ca)*2);
  STDERR_BASE_LINE("sf=%g,ni=%d,ca=%g deg,log2ca_x2=%d\n",sf,ni,ca,log2ca_x2);
  ca=exp2(log2ca_x2/2)*(log2ca_x2&1?1.5:1);
  if(mesh!=0){
    STDERR_BASE_LINE("ca=%g\n",ca);
    ca*=.5;
    }
  
  const double cd2=square(ca*d2r*sf);
  STDERR_BASE_LINE("sf=%g,ca=%g deg,cd=%g\n",sf,ca,sqrt(cd2));
  
  const int cny=180./ca,cnx=360./ca,cn=cnx*cny;
  int ci,cj,cm;/*< clustering indexes */
  
  static vector3_t *cu=0;
  static T *cmass=0;/*< mass of the clusters */
  static vector<int> *cIs=0;
  bool doCIs=false;
  
  if(cu==0 or cmass==0 or cIs==0){
    
    if(cu!=0 or cmass!=0 or cIs!=0) TRAP_ERR_EXIT(ENOEXEC,"programming error: %p!=0 or %p!=0 or %p!=0\n",cu,cmass,cIs);
    
    cu=new vector3_t[cn];
    
    for(cm=0;cm<cn;cm++){
      cj=cm/cnx;
      ci=cm-cj*cnx;
      
      const double lat=(cj+.5)*ca-90.;
      const double lon=(ci+.5)*ca;
      cLat=cos(lat*d2r);
      sLat=sin(lat*d2r);
      cLon=cos(lon*d2r);
      sLon=sin(lon*d2r);
      
      vector3_t *cucm=&cu[cm];
      cucm->init(cLat*cLon,cLat*sLon,sLat);
      }
    
    cmass=new T[cn];
    cIs=new vector<int>[cn];
    doCIs=true;
    }
  
  aset(cmass,cn,(T)0.);
  
  for(ii=0;ii<ni;ii++){
    double loii=loi[ii];
    degree_recale(&loii,180.);
    
    ci=floor(loii/ca);
    cj=floor((lai[ii]+90.)/ca);
    cm=cj*cnx+ci;
    
    if(cm<0 or cn<=cm){
      if(mesh==0 and mass[ii]==0.)
        continue;
      TRAP_ERR_EXIT(ENOEXEC,"programming error: cm=%d*%d+%d=%d not in [0;%d[\n",cj,cnx,ci,cm,cn);
      }
    if(doCIs)
      cIs[cm].push_back(ii);
    cmass[cm]+=mass[ii];
    }
  
#if CLUSTERING_ENABLED && 0
  if(debugPath!=0){
    /* using another path as mesh and grid are likely to interfere */
    string path=replace(debugPath,".nc","-clusters.nc");
    unlink(path.c_str());
    
    int status;
    grid_t cgrid;
    const double hca=ca*.5;
    map_set2Dgrid(&cgrid,hca,-90+hca,360-hca,90-hca,ca);
    map_completegridaxis(&cgrid,1);
    poc_var_t avar,gvar;
    status=poc_save_grid(path,&avar,cgrid,0,1);
    avar.init("mass_a",NC_DOUBLE,"","kg");
    (gvar=avar).init("mass_G",NC_DOUBLE,"","degrees");
    status=poc_put_cvara(path,avar,gvar,0,cmass);
    cgrid.free();
    
    status=poc_def_att(path,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
    status=poc_def_att(path,poc_att_t("history",cmd));
    }
#endif
  STDOUT_BASE_LINE("Made %dx%d=%d clusters in %g s.\n",cnx,cny,cn,difftime(before));
#endif
  
/*------------------------------------------------------------------------------
  calculate loading effect */
  /// \bug Damien Allain 2012-07-24 : assuming constant gravity when it varies by +/-0.3%. See gravitational_constant_01().
  
// #warning comment this line and the right ones below
  range_t<int> ioR
    (0,no);
    //(no-no/60,no);
    //(0,no/60);
    //(58e4,61e4);
    //(35e4,45e4);
  ioR.narrow(0,no);
  const int ioRd=ioR.d();
  if(0<ioR.min or ioR.max<no){
    STDOUT_BASE_LINE("*** Only %d(%g%%) of %d points, from %d to %d ***\n",ioRd,100.*ioRd/no,no,ioR.min,ioR.max);
    ioR.narrow(0,no-1);
    }
  else
    STDOUT_BASE_LINE("%d points\n",no);
  range_t<double> loR
    (-400,400);
    //(-20,10);
    //(60,90);
  
  STDERR_BASE_LINE("PID: %d\n",getpid());
  gettimeofday(&before);
  
  if(mesh!=0){
    int *P1_ordinals=new int[no];
    
    for(io=0;io<no;io++)
      P1_ordinals[io]=io;
    
    poc_put_mesh_vara(outputPath,*mesh,LGP1,"P1_ordinals",-1,P1_ordinals,"",0);
    
    delete[] P1_ordinals;
    }
  
  int nprocs=initialize_OPENMP();
  
  #pragma omp parallel num_threads(nprocs)
  {
  int cm,ii,io,it;
  
  #pragma omp for schedule(dynamic,1)
  for(io=0;io<no;io++){
    
    if(io==0){
      STDOUT_BASE_LINE("omp_get_num_threads()=%d . omp_get_dynamic()=%d .\n",omp_get_num_threads(),omp_get_dynamic());
      fflush(stdout);
      }
    
    //loop and header take 0.02s
    const vector3_t uoi=uo[io];
    const double ek=1/hypot(uoi.x,uoi.y);
    const vector3_t
      eoi(-uoi.y*ek,uoi.x*ek,0.),
      noi=uoi&eoi;
#if FAST_GRAD
    const double
      fastGradD=.18,
      fastGradD2=square(fastGradD);
    const vector3_t
      fastGradV(fastGradD,fastGradD,fastGradD),
      uMin=uoi-fastGradV,
      uMax=uoi+fastGradV;
#endif
    double d2;/*< square of distance on unit sphere */
    T loadio=0.,potio=0.,lsaio=0.;
    T eastio=0.,nortio=0.;
    
/*------------------------------------------------------------------------------
    assertion and tests */
#if 1
    {//0.22s, REGARDLESS OF THE NUMBER OF NODES
    time_t now;
    if((now=timeIsOld()) or io==ioR.max){
      const int ioD=io-ioR.min;
      double percent=100.*ioD/ioRd;
      fprintf(stderr,"%s%s%d/%d;%.03g%% in %.03gs",cr,el,ioD,ioRd,percent,difftime(before));
#if CLUSTERING_ENABLED && 0
      const int
        period=0x20,
        periodMask=period-1;
      static double oldPercent=0.;
      if(true and mesh!=NULL){
        /* progressive save */
        if(not (now&periodMask) and 1>oldPercent-percent){
          oldPercent=percent;
          int status,varStatus;
          fprintf(stderr," (Saving to %s ...",outputPath);
          struct timeval saveStart;
          gettimeofday(&saveStart);
          varStatus=poc_inq_var(outputPath,"LSA_a",0,-1);
          if(not lsaOnly){
            status=poc_put_mesh_cvara(outputPath,*mesh,LGP1,"UP_a","UP_G",0,displacement,"m",0);
            if(deformation!=0)
              status=poc_put_mesh_cvara(outputPath,*mesh,LGP1,"UP_CF_a","UP_CF_G",0,deformation,"m",0);
#if HORIZONTAL_COLUMN >0
            status=poc_put_mesh_cvara(outputPath,*mesh,LGP1,"EAST_a","EAST_G",0,east,"m",0);
            status=poc_put_mesh_cvara(outputPath,*mesh,LGP1,"NORTH_a","NORTH_G",0,north,"m",0);
#endif
            status=poc_put_mesh_cvara(outputPath,*mesh,LGP1,"potential_a","potential_G",0,pot,"m",0);
            }
          status=poc_put_mesh_cvara(outputPath,*mesh,LGP1,"LSA_a","LSA_G",0,lsa,"m",0);
          if(varStatus!=0){
            status=poc_def_att(outputPath,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
            status=poc_def_att(outputPath,poc_att_t("history",cmd));
            }
          fprintf(stderr,"%g s)",difftime(saveStart));
          }
        else
          fprintf(stderr," (Saving in %ds)",period-now&periodMask);
        }
#endif
      }
    if(io==ioR.max) STDERR_BASE_LINE("*** TESTING ***");
    }
#endif
    if(not ioR.has(io) or not loR.has(degree_recale(loo[io],0.)) ){//<0s on full set ?!, 0.02s on 10% of full set.
      
      if(not lsaOnly){
        
        if(deformation!=0)
          deformation[io]=(T)NAN;
        
        displacement[io]=(T)NAN;
#if HORIZONTAL_COLUMN > 0
        east[io]=(T)NAN;
        north[io]=(T)NAN;
#endif
        pot[io]=(T)NAN;
        }
      
      lsa[io]=(T)NAN;
      
      continue;
      }
    
    /* time budgets are for executables compiled without -O3 */
    
#if CLUSTERING_ENABLED
/*------------------------------------------------------------------------------
    for each cluster */
    for(cm=0;cm<cn;cm++){/* 5.4s */
      
//       #warning comment this line and the one below
//       if(cm*2!=cn)continue; /* test only one cluster */
      
      const vector3_t *cucm=&cu[cm];/* 4.6-5.4=-0.8s (!) */
      d2=square(*cucm-uoi);/* 45s */
      
      if(d2>cd2){/* 4.6s */
        const T *cmasscm=&cmass[cm];/* 3.5s */
        if(*cmasscm==0.) continue;/* 9.3s */
        it=tableInterpolation*(maxFlog2r2-flog2(d2));/* 45s */
        lsa_increment(table,it,*cmasscm,eoi,noi,*cucm,&loadio,&eastio,&nortio,&potio,&lsaio);
        continue;
        } /* 156s or 310s */
      
      const vector<int> &cIscm=cIs[cm];
      const size_t cIscmL=cIscm.size();
      int cIscmi;
      
/*------------------------------------------------------------------------------
      for each input element */
      for(cIscmi=0;cIscmi<cIscmL;cIscmi++)
#else
      for(ii=0;ii<ni;ii++)//2.2min
#endif
        {
#if CLUSTERING_ENABLED
//         #warning comment this line and the one below
//         if(grid!=0 and cIscmi*2!=cIscmL)continue; /* test only one grid element */
        
        ii=cIscm[cIscmi];
#endif
        
        if(mass[ii]==0.)
          continue;
        
        const vector3_t *uii=&ui[ii];
#if FAST_GRAD
        if(
            uii->x<uMin.x or uMax.x<uii->x or
            uii->y<uMin.y or uMax.y<uii->y or
            uii->z<uMin.z or uMax.z<uii->z)
          continue;
#endif
        d2=square(*uii-uoi);//2min
#if FAST_GRAD
        if(d2>fastGradD2)
          continue;
#endif
        
/*------------------------------------------------------------------------------
        if the mesh is given */
        if(mesh!=NULL && d2<maxFracD2){
          
          /* check distances of each corner */
          range_t<double> d2R;
          int i,ic;
          for(i=0;i<3;i++){
            ic=mesh->triangles[ii].vertex[i];
            d2R << square(uo[ic]-uoi);
            }
/*------------------------------------------------------------------------------
          if the ouput node is one of the corners of the input triangle */
          if(d2R.min==0.){
            int Is[2],Ii=0,&i1=Is[0],&i2=Is[1];
            i1=-1;i2=-1;
            
            /* get other corners */
            for(i=0;i<3;i++){
              if(mesh->triangles[ii].vertex[i]==io)
                continue;
              Is[Ii]=mesh->triangles[ii].vertex[i];
              Ii++;
              }
            if(i1<0 or i2<0)
              TRAP_ERR_EXIT(ENOEXEC,"programming error: %d and %d\n",i1,i2);
            
            /* integrate */
            const int N=8;
            vector3_t u=uo[i1];
            const double k=1./N;
            const vector3_t du=(uo[i2]-u)*k;
            const T massik=mass[ii]*k;
            
            u+=du*.5;
            
            for(i=0;i<N;i++,u+=du){
              d2=square(uoi-u);
              it=tableInterpolation*(maxFlog2r2-flog2(d2));
              if(it>=tableLength)
                TRAP_ERR_EXIT(ENOEXEC,"programming error: [(%g;%g)->(%g;%g)[%d/%d]]=%g<%g\n",
                  loi[ii],lai[ii],loo[io],lao[io],i,N,
                  sqrt(d2),sqrt(shortestEgdeL2));
              lsa_increment(CADTable,it,massik,eoi,noi,u,&loadio,&eastio,&nortio,&potio,&lsaio);
              }
            
            continue;
            }
  /*------------------------------------------------------------------------------
          if the ouput node is close to the input element, relatively to the size of the latter */
          if(d2R.min<0.8*d2R.max){
            /* integrate */
            int j;
            double s;
            
            int N=ceil((d2R.max/d2R.min-1)*15);
            /* N can be extremely large (>1e5) for a point e.g. on the other
            side of a breakwater and so cause an almost infinite loop */
// #warning comment this line and put value below back to 100
            N=min(N,100);
            const double k=1./N;
            const T massik2=mass[ii]*square(k);
            
            const int *Is=mesh->triangles[ii].vertex;
            const vector3_t
              *u0=&uo[Is[0]],
              *u1=&uo[Is[1]],
              *u2=&uo[Is[2]],
              dui=(*u1-*u0)*k,
              duj=(*u2-*u0)*k,
              dup=(dui+duj)*(1./3)-dui*.5;
            vector3_t u,
              uj=*u0+dui*.5;
            
            for(j=0;j<N;j++,uj+=duj){
              for(s=1.;s>=-1.;s-=2.){
                u=uj+dup*s;
                for(i=0;i<N-j;i++,u+=dui){
                  d2=square(uoi-u);
                  it=tableInterpolation*(maxFlog2r2-flog2(d2));
                  lsa_increment(table,it,massik2,eoi,noi,u,&loadio,&eastio,&nortio,&potio,&lsaio);
                  }
                
                if(j==0)
                  break;
                }
              }
            
            continue;
            }
          }//2.5min
/*------------------------------------------------------------------------------
        if the grid is given and the node is close */
        if(grid!=NULL and d2<maxFracD2){
          
          const int
            j0=ii/grid->nx,
            i0=ii-j0*grid->nx;
          int di=-1;
          if(i0==0) di+=grid->nx-(grid->nx&1);
          
          const vector3_t
            *u0=&uo[ii],
            *u1=&uo[ii+di],
            *u2=&uo[ii-grid->nx];
          
          int i,j;
          
          const bool
            lineIntegrate=(ii==io) or (i0==0 and (grid->nx&1) and io-ii==grid->nx-1);
          
          int N;
          if(lineIntegrate)
            N=8;
          else
            N=sqrt(ceil(maxFracD2/d2));
          const double k=1./N;
          const vector3_t
            dui=(*u0-*u1)*k,
            duj=(*u0-*u2)*k;
          vector3_t u;
          
          /* if the input node is on the output node */
          if(lineIntegrate){
            
            /* line-integrate */
            const T massik=mass[ii]*(0.25*k);
            
            vector3_t du;
            int sign;
            
            for(j=0;j<4;j++){
              
              sign=(j&2)?-1:1;
              du=sign*((j&1)?duj:dui);
              u=*u0-((dui+duj)*(sign*N)-du)*.5;
              
              for(i=0;i<N;i++,u+=du){
                d2=square(uoi-u);
                it=tableInterpolation*(maxFlog2r2-flog2(d2));
                if(it>=tableLength) TRAP_ERR_EXIT(ENOEXEC,"programming error: %g<%g\n",sqrt(d2),sqrt(shortestEgdeL2));
                lsa_increment(CADTable,it,massik,eoi,noi,u,&loadio,&eastio,&nortio,&potio,&lsaio);
                }
              }
            
            continue;
            }
          
          /* otherwise surface integrate */
          const T massik2=mass[ii]*square(k);
          vector3_t uj;
          
          uj=*u0-(dui+duj)*((N-1)*0.5);
          
          for(j=0;j<N;j++,uj+=duj){
            u=uj;
            for(i=0;i<N;i++,u+=dui){
              d2=square(uoi-u);
              it=tableInterpolation*(maxFlog2r2-flog2(d2));
              lsa_increment(table,it,massik2,eoi,noi,u,&loadio,&eastio,&nortio,&potio,&lsaio);
              }
            }
          
          continue;
          }
        
/*------------------------------------------------------------------------------
        if the output node is far enough from the input element, relatively to the size of the latter */
        it=tableInterpolation*(maxFlog2r2-flog2(d2));//5min
        if(it>=tableLength)
          #pragma omp critical(human)
          TRAP_ERR_EXIT(ENOEXEC,"programming error for (%g;%g)->(%g;%g): %g<%g\n",
            loi[ii],lai[ii],loo[io],lao[io],sqrt(d2),sqrt(shortestEgdeL2));
        const T *massi=&mass[ii];
        lsa_increment(table,it,*massi,eoi,noi,*uii,&loadio,&eastio,&nortio,&potio,&lsaio);
        
        }
      
#if CLUSTERING_ENABLED
      }
#endif
    
#if POT_UP_TO_LSA == 1
    lsaio=potio-loadio;
#elif POT_UP_TO_LSA == -1
#else
#error POT_UP_TO_LSA not properly defined
#endif
    
    if(not lsaOnly){
#if GREEN_MODE == GREEN_DEFORMATION
      if(deformation!=0)
        deformation[io]=(T)loadio;
      displacement[io]=(T)(loadio-momentumDisplacement*uoi);
#if HORIZONTAL_COLUMN > 0
#error not coded yet
#endif
#elif GREEN_MODE == GREEN_DISPLACEMENT
      if(deformation!=0)
        deformation[io]=(T)(loadio+momentumDisplacement*uoi);
      displacement[io]=(T)loadio;
#if HORIZONTAL_COLUMN > 0
      east[io]=(T)eastio;
      north[io]=(T)nortio;
#endif
#else
#error GREEN_MODE not properly defined
#endif
      
#if POT_UP_TO_LSA == -1
      potio=lsaio+loadio;
#endif
      pot[io]=(T)potio;
      }
    lsa[io]=(T)lsaio;
    }
  
  }/* EO pragma omp parallel */
  double cct=difftime(before);//core computation time
  fprintf(stderr,"\n");
  STDOUT_BASE_LINE("Took %g s.\n",cct);
  
  return cct;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void print_help(const char *prog_name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/** \brief prints help of this programme.

\param prog_name name of the program to be printed by this help function : argv[0]
*/
/*----------------------------------------------------------------------------*/
{
/** The code of the body of this function is :
    \code /**/ // COMPILED CODE BELOW !!!
  printf("\n"
    "NAME AND VERSION\n"
    "   %s version " VERSION " " REVISION "\n", prog_name);
  printf("\n"
    "USE\n"
    "  %s [OPTIONS]\n",prog_name);
  printf("\n"
    "DESCRIPTION\n"
/* USE FIND AND REPLACE TO EDIT THE LINE BELOW, SO THAT THE FILE HEAD IS ALSO UPDATED */
    "  Computes the tidal loading.\n"
    "  Outputs the radial deformation.\n"
    "\n"
    "OPTIONS :\n"
    "  -z  followed by either :\n"
    "        regular grid name\n"
    "        or regular global grid resolution in degrees\n"
    "        or unstructured netcdf file.\n"
    "  -m  followed by continental unstructured mesh on which to interpolate the unstructured netcdf file. If this last one is not given, test algorithm with 1m uniform input.\n"
    "  -r  followed by regular netcdf elevations.\n"
    "  -v  obsolete. Use -e instead\n"
    "  -e  followed by 2 variable names for amplitude and phase of elevation. Default is Ha Hg\n"
    "  -u  followed by 4 variable names for amplitude and phase of zonal and meridional speeds. If not given, the angular momentum is not calculated\n"
    "  -b  followed by topography file. Usefull for masking.\n"
    "  -o  followed by output file name. If omitted, stop after mass budget.\n"
    "  -D  followed by debug output file name, where re-discretised elevations, etc... will be saved\n"
    "  --lsa-only : only produce LSA to save time. Please note altimetry actually needs displacements!\n"
    "\n"
    "EXAMPLES\n"
    "loading -z /data1/FES2014/M2.EnOI.nc -e analysis_Ha analysis_Hg -m FES2014-continents.nei -o M2.EnOI-LSA-continents.nc\n"
    "loading -r /data1/FES2014/soa/M2.FES2014b.nc -e elevation_a elevation_G -o M2.FES2014b-LSA.nc\n"
    );
  print_OPENMP_help(prog_name); /** \endcode */
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_get_triangle_lon_lat(const mesh_t & mesh,int t,double *lon,double *lat)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int c,n;/*< indexes : corners, vertices */
  double y=0.,x=0.,x0=NAN;
  
  const triangle_t *triangle=&mesh.triangles[t];
  const vertex_t *vertex;
  
  for(c=0;c<3;c++){
    n=triangle->vertex[c];
    vertex=&mesh.vertices[n];
    
    y+=vertex->lat;
    
    if(c==0){
      x0=vertex->lon;
      x+=x0;
      }
    else
      x+=degree_recale(vertex->lon,x0);
    }
  
  y/=3.;
  x/=3.;
  
  *lon=x;
  *lat=y;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T *toP0(const char *descStr,mesh_t *mesh,int discretisation,T *data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const char *discname;
  const discretisation_t *descriptor=0;
  
  discretisation_init(mesh,discretisation,-1);
  descriptor=get_descriptor_address(*mesh,discretisation);
  discname=discretisation_name(discretisation);
  STDOUT_BASE_LINE("%s input : %d %s nodes\n",descStr,descriptor->nnodes,discname);
  
  int c0,cm,cd;
  
  switch(discretisation){
  case LGP1:
  case NCP1:
  case DGP1:
  case DNP1:
    c0=0;
    cm=3;
    cd=1;
    break;
  case LGP2:
    c0=1;
    cm=6;
    cd=2;
    break;
  default:
    TRAP_ERR_RETURN((T *)0,1,"not coded for %s discretisation\n",discname);
    }
  
  T *P0=new T[mesh->ntriangles];
  
  T interpolated;
  double w;
  int t,c,n;//indexes : triangles, corners, nodes
  
  for(t=0;t<mesh->ntriangles;t++){
    w=0;
    interpolated=0.;
    for(c=c0;c<cm;c+=cd){
      n=descriptor->NIbE[t][c];
      interpolated+=data[n];
      w++;
      }
    if(w>1)
      interpolated/=w;
    
    if( isnan(interpolated) ){
      double lon,lat;
      
      fe_get_triangle_lon_lat(*mesh,t,&lon,&lat);
      
      TRAP_ERR_EXIT(-1,
        "nan for %d/%d at (%g;%g): %g\n",
        t,mesh->ntriangles,lon,lat,interpolated);
      }
    P0[t]=interpolated;
    }
  
  return P0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void fe_get_true_areas(const mesh_t & mesh, double *lon, double *lat, double *area, double *perimeter_area)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int t,n,c;
  int areaErrorC=0;
  
  for(t=0;t<mesh.ntriangles;t++){
    vector3_t v;
    const triangle_t *triangle=&mesh.triangles[t];
    const vertex_t *vertex,*corners[3];
    
    /* coordinates */
    v.init(0.,0.,0.);
    
    for(c=0;c<3;c++){
      n=triangle->vertex[c];
      vertex=&mesh.vertices[n];
      
      v+=math_polar2cartesian(vertex->lon*d2r,vertex->lat*d2r);
      
      corners[c]=vertex;
      }
    
    double lont,latt;
    math_cartesian2polar(v,&lont,&latt);
    lont*=r2d;
    latt*=r2d;
    lon[t]=lont;
    lat[t]=latt;
    
    /* proper area */
    double areat,cornerAngle[3],areaError;
    
    /* Demonstration of the solid angle formula:
    Given a triangle ABC,
    if you fill up the angles A, B and C with their opposites,
    you will have filled-up:
    + 3 times the triangle and its opposite
    + once the remaining sphere
    So: 4*pi+(3-1)*2*area=4*(A+B+C)
    or: 4*area=4*(A+B+C-pi)
    */
    
    for(c=0;c<3;c++){
      const int
        c1=(c+1)%3,
        c2=(c+2)%3;
      
      cornerAngle[c]=fabs(geo_angle_radian(corners[c]->lon,corners[c]->lat,
        corners[c1]->lon,corners[c1]->lat,
        corners[c2]->lon,corners[c2]->lat));
      }
    
    areat=(cornerAngle[0]+cornerAngle[1]+cornerAngle[2])-M_PI;
    areat*=MeanEarthRadius2;
    
    areaError=fabs(triangle->TrueArea/areat-1);
    if(areaError>.1){
      areaErrorC++;
      if(areaErrorC==1)
        STDERR_BASE_LINE("WRONGEST AREAS:\n");
      STDERR_BASE_LINE("[%d](%g;%g): %g <-> %g (%d)\n",t,lont,latt,triangle->TrueArea,areat,areaErrorC);
      }
    
    area[t]=areat;
    
    /* perimeter / area ratio */
    if(perimeter_area==0) continue;
    
    edge_t *edge;
    double perimeter=0;
    
    for(c=0;c<3;c++){
      n=triangle->edges[c];
      edge=&mesh.edges[n];
      if(edge->L==0.){
        n=edge->extremity[0];
        vertex=&mesh.vertices[n];
        
        n=edge->extremity[1];
        
        edge->L=geo_distance(vertex->lon,vertex->lat,mesh.vertices[n].lon,mesh.vertices[n].lat);
        }
      perimeter+=edge->L;
      }
    
    perimeter_area[t]=perimeter/triangle->TrueArea;
    }
  
}


#define T double
#define NC_T NC_DOUBLE
#define NC_FILL_T NC_FILL_DOUBLE
#define TS "double"
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loading_UG()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// tidal loading with an unstructured grid
/*----------------------------------------------------------------------------*/
{
STDERR_BASE_LINE("using " TS " for saving results\n");
  
  int t,n;//indexes : triangles, corners, nodes
  double *lon,*lat;//triangle coordinates
  double *outputLon,*outputLat;//node coordinates

  int status;
  mesh_t *inputMesh=0,*outputMesh=0;
  double *area,*perimeter_area=0;

/*------------------------------------------------------------------------------
  meshes and grid input */
  
  if(topoPath!=NULL)STDERR_BASE_LINE("ignoring topography from %s\n",topoPath);
  
  if(meshPath!=0){
    outputMesh=new mesh_t[1];
    status=fe_readmesh(meshPath,MESH_FILE_FORMAT_UNKNOWN,outputMesh);
    status=fe_list(outputMesh);//elements
    status=fe_edgetable(outputMesh,0,0);
    }
  
  bool isComplex= (elevPhaseName!=0 and *elevPhaseName!='\0');
  int frame;
  
  if(isComplex)
    frame=-1;
  else
    frame=0;
  
  poc_var_t ampVar;
  int discretisation;
  
  if(zone!=0){
    
    status=poc_inq_var(zone,elevAmplitudeName,&ampVar);
    if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_inq_var(\"%s\",\"%s\",) error", zone,elevAmplitudeName);
    
    inputMesh=new mesh_t[1];
    status=poc_get_mesh(zone,ampVar,inputMesh,1,frame,&discretisation);
    if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_get_mesh error with %s",zone);
    }
  
  grid_t listGrid,grid;
  list_t list;
  string SP_ampName;
  
  if(outputMesh==0 and inputMesh!=0){
    STDOUT_BASE_LINE("Using input mesh for output\n");
    outputMesh=inputMesh;
    }
  else if(outputMesh!=0 and inputMesh==0){
    }
  else if(outputMesh!=0 and inputMesh!=0){
    fe_Allocate_and_CreateList(*inputMesh,&listGrid,&list);
    
    if(inputPath!=0){
      
      if(isComplex){
        SP_ampName="sp_a";
        }
      else{
        SP_ampName="sp";
        }
      
      if( access(inputPath,R_OK)==0 ){
        status=poc_inq_var(inputPath,SP_ampName,&ampVar);
        if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_inq_var(\"%s\",\""+SP_ampName+"\",) error", inputPath);
        
        status=poc_get_grid(inputPath,ampVar,&grid);
#warning hard coded
        if(status!=NC_NOERR)
          status=poc_get_grid("/data1/.data3/ECMWF/EI_AN+FC/EI_AN+FC_200001.EC",ampVar,&grid);
        if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_get_grid(\"%s\",(\""+ampVar.name+"\"),) error",inputPath);
        }
      
      }
    
    }
  else TRAP_ERR_EXIT(ENOEXEC,"INPUT IS %s:%p, OUTPUT IS %s:%p\n",zone,inputMesh,meshPath,outputMesh);
  
  STDOUT_BASE_LINE("OUTPUT: %d triangles and %d points... ",outputMesh->ntriangles,outputMesh->nvtxs);
  fflush(stdout);
  
  /* centers of mass of input triangles */
  area=new double[outputMesh->ntriangles];
  if(debugPath!=0)
    perimeter_area=new double[outputMesh->ntriangles];
  lat=new double[outputMesh->ntriangles];
  lon=new double[outputMesh->ntriangles];
  
  printf("allocated\n");fflush(stdout);
  
  fe_get_true_areas(*outputMesh,lon,lat,area,perimeter_area);
  
  outputLat=new double[outputMesh->nvtxs];
  outputLon=new double[outputMesh->nvtxs];
  for(n=0;n<outputMesh->nvtxs;n++){
    const vertex_t *vertex=&outputMesh->vertices[n];
    outputLat[n]=vertex->lat;
    outputLon[n]=vertex->lon;
    }
  
/*------------------------------------------------------------------------------
  field input */
  
  if(isComplex){
    
    complex<double> *elevation=0;
    complex<double> *masses;
    complex<T> *deformation=0,*displacement=0,*east=0,*north=0,*pot=0,*lsa;
    
    if(zone!=0){
      poc_cdata_t input;
      
      status=input.init(zone,elevAmplitudeName,elevPhaseName,"T",indexes_t()<<frame);
      if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_cdata_t::init(\"%s\",\"%s\",\"%s\",...) error", zone,elevAmplitudeName,elevPhaseName);
      
      elevation=toP0("elevation",inputMesh,discretisation,input.data);
      
      if( uAmpName!=0 and uPhaName!=0 and vAmpName!=0 and vPhaName!=0 ){
        status=input.init(zone,uAmpName,uPhaName,"T",indexes_t()<<frame);
        if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_cdata_t::init(\"%s\",\"%s\",\"%s\",...) error",zone,uAmpName,uPhaName);
        
        discretisation=poc_get_discretisation(input.info,inputMesh);
        uTide=toP0("u",inputMesh,discretisation,input.data);
        
        status=input.init(zone,vAmpName,vPhaName,"T",indexes_t()<<frame);
        if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_cdata_t::init(\"%s\",\"%s\",\"%s\",...) error",zone,vAmpName,vPhaName);
        
        discretisation=poc_get_discretisation(input.info,inputMesh);
        vTide=toP0("v",inputMesh,discretisation,input.data);
        
        poc_data_t<double> depthInput;
        string depthFile=zone;
        
        if( strrncmp(zone,".EnOI.nc")==0 ){
          replace(&depthFile,"EnOI","spectral");
          STDOUT_BASE_LINE("*** PATCHING BATHYMETRY INPUT FILE TO "+depthFile+" ***\n");
          }
        
        status=depthInput.init(depthFile,"bathymetry");
        if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_data_t::init(\""+depthFile+"\",\"bathymetry\") error");
        
        discretisation=poc_get_discretisation(depthInput.info,inputMesh);
        depth=toP0("depth",inputMesh,discretisation,depthInput.data);
        
        depthInput.destroy_data();
        }
      
      input.destroy_data();
      
      if(frame<=0 and debugPath!=0 and inputMesh!=0){
        /* using another path as the input mesh may be different from the output mesh */
        string path=replace(debugPath,".nc","-0.nc");
        unlink(path.c_str());
        
        int *P0_ordinals=new int[inputMesh->ntriangles];
        for(t=0;t<inputMesh->ntriangles;t++)
          P0_ordinals[t]=t;
        
        status=poc_put_mesh_vara(path.c_str(),*inputMesh,LGP0,"P0_ordinals",-1,P0_ordinals,"",0);
        NC_CHKERR_BASE_LINE(status,"poc_put_mesh_vara(\""+path+"\",,,\"P0_ordinals\",...) error");
        
        delete[] P0_ordinals;
        
        status=poc_put_mesh_cvara(path.c_str(),*inputMesh,LGP0,"a_eta_P0","G_eta_P0",-1,elevation,"m",1);
        NC_CHKERR_BASE_LINE(status,"poc_put_mesh_cvara(\""+path+"\",,,\"a_eta_P0\",\"G_eta_P0\",...) error");
        
        status=poc_def_att(path,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
        status=poc_def_att(path,poc_att_t("history",cmd));
        }
      }
    
    if(outputMesh!=0 and inputMesh==0){
      if(frame>0) TRAP_ERR_EXIT(-1,"UNIFORM 1m INPUT not relevant for time series\n");
      STDOUT_BASE_LINE("FAKING UNIFORM 1m INPUT\n");
      elevation=aset(outputMesh->ntriangles,complex<double>(1.,0.));
      }
    else if(outputMesh!=0 and inputMesh!=0 and outputMesh!=inputMesh){
      STDOUT_BASE_LINE("reinterpolating\n");
      int last=-1;
      int64_t gridLast=-1;
      
      complex<double> *oldElevation=elevation;
      if(uTide!=0 or vTide!=0){
        STDERR_BASE_LINE("WARNING: SPEED REINTERPOLATION NOT CODED YET: CANCELING ANGULAR MOMENTUM.");
        deletep(&uTide);
        deletep(&vTide);
        deletep(&depth);
        }
      elevation=new complex<double>[outputMesh->ntriangles];
      const complex<double> mask=0.;/* MUST BE 0. */
      
      poc_cdata_t surfacePressure;
      
      if(inputPath!=0 and SP_ampName!=""){
        const string SP_phaName="sp_G";
        
        status=surfacePressure.init(inputPath,SP_ampName,SP_phaName,"T",indexes_t()<<frame);
        if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_cdata_t::init(\"%s\",\""+SP_ampName+"\",\""+SP_phaName+"\") error",inputPath);
        }
      
      for(t=0;t<outputMesh->ntriangles;t++){
        double x,y; /*< interpolation coordinates */
        complex<double> interpolated;
        int verbose=0;
        
        fe_get_triangle_lon_lat(*outputMesh,t,&x,&y);
        
        //verbose= t==397747;
        if(verbose)
          STDERR_BASE_LINE("[%d](%g;%g):",t,x,y);
        
        n=fe_detect(*inputMesh,inputMesh->triangles,list,listGrid,x,y,last);
        status=fe_interpolate2D(*inputMesh,LGP0,oldElevation,mask,x,y,n,&interpolated);
        if(verbose) fprintf(stderr,"[%d] %d %g,%g°\n",n,status,abs(interpolated),-arg(interpolated)*r2d);
        
        if(isnan(interpolated)) TRAP_ERR_EXIT(ENOEXEC,"NAN found at [%d](%g,%g)\n",t,x,y);
        
        if(surfacePressure.data!=0){
          complex<double> z;
          updatemin(&y,grid.ymax);
          updatemax(&y,grid.ymin);
          index_interpolation(grid,x,y,&gridLast,surfacePressure.data,mask,&z);
          if(isnan(z)) TRAP_ERR_EXIT(ENOEXEC,"NAN found at [%d](%g,%g)\n",t,x,y);
          z*=1./(1025*9.81);
          interpolated+=z;
          }
        
        elevation[t]=interpolated;
        }
      
      delete[]oldElevation;
      }
    else if(outputMesh==inputMesh){
      }
    else TRAP_ERR_EXIT(ENOEXEC,"programming error: should have exited before : inputMesh=%p; outputMesh=%p;\n",inputMesh,outputMesh);
    
    if(inputMesh!=0 and inputMesh!=outputMesh){
      inputMesh->destroy();
      deletep(&inputMesh);
      }
    
    /** \bug 2012-05-22 Damien Allain : using a constant volumic mass */
    const double rho_water=1025.;/*< volumic mass of sea water */
    double waterColumMass;
    
    if(uTide!=0 and vTide!=0 and depth!=0){
      for(t=0;t<outputMesh->ntriangles;t++){
      /* local momentum */
        waterColumMass=area[t]*rho_water*depth[t];
        uTide[t]*=waterColumMass;
        vTide[t]*=waterColumMass;
        }
      }
    
    deletep(&depth);
    
    if(frame<=0 and debugPath!=0){
/*------------------------------------------------------------------------------
      save output mesh */
      unlink(debugPath);
      
      int *P0_ordinals=new int[outputMesh->ntriangles];
      for(t=0;t<outputMesh->ntriangles;t++)
        P0_ordinals[t]=t;
      
      status=poc_put_mesh_vara(debugPath,*outputMesh,LGP0,"P0_ordinals",-1,P0_ordinals,"",0);
      NC_CHKERR_BASE_LINE(status,"poc_put_mesh_vara(\"%s\",,,\"P0_ordinals\",...) error",debugPath);
      
      delete[] P0_ordinals;
      
      status=poc_put_mesh_cvara(debugPath,*outputMesh,LGP0,"a_eta_P0","G_eta_P0",-1,elevation,"m",1);
      status=poc_put_mesh_vara(debugPath,*outputMesh,LGP0,"area",-1,area,"m2",1);
      status=poc_put_mesh_vara(debugPath,*outputMesh,LGP0,"perimeter_area",-1,perimeter_area,"m",1);
      status=poc_def_att(debugPath,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
      status=poc_def_att(debugPath,poc_att_t("history",cmd));
      }
    
/*------------------------------------------------------------------------------
    complex masses */
    masses=massFromAreaAndElevation<complex<T> >(outputMesh->ntriangles,area,elevation);
    delete[]elevation;
    deletep(&area);
    
/*------------------------------------------------------------------------------
    call loading() */
    
    //deformation=new complex<T>[outputMesh->nvtxs];
    if(not lsaOnly){
      displacement=new complex<T>[outputMesh->nvtxs];
  #if HORIZONTAL_COLUMN >0
      east=new complex<T>[outputMesh->nvtxs];
      north=new complex<T>[outputMesh->nvtxs];
  #endif
      pot=new complex<T>[outputMesh->nvtxs];
      }
    lsa=new complex<T>[outputMesh->nvtxs];
    
    if(frame<=0)
      unlink(outputPath);
    
    loading<cvector3_t,zmatrix3x3_t>(outputMesh->ntriangles,lat,lon,masses,outputMesh->nvtxs,outputLat,outputLon,deformation,displacement,east,north,pot,lsa,outputMesh);
    
/*------------------------------------------------------------------------------
    save results */
    if(not lsaOnly){
      status=poc_put_mesh_cvara(outputPath,*outputMesh,LGP1,"UP_a","UP_G",0,displacement,"m",1);
      if(deformation!=0)
        status=poc_put_mesh_cvara(outputPath,*outputMesh,LGP1,"UP_CF_a","UP_CF_G",0,deformation,"m",1);
  #if HORIZONTAL_COLUMN >0
      status=poc_put_mesh_cvara(outputPath,*outputMesh,LGP1,"EAST_a","EAST_G",0,east,"m",1);
      status=poc_put_mesh_cvara(outputPath,*outputMesh,LGP1,"NORTH_a","NORTH_G",0,north,"m",1);
  #endif
      status=poc_put_mesh_cvara(outputPath,*outputMesh,LGP1,"potential_a","potential_G",0,pot,"m",1);
      }
    status=poc_put_mesh_cvara(outputPath,*outputMesh,LGP1,"LSA_a","LSA_G",0,lsa,"m",1);
    
    if(frame<=0){
      status=poc_def_att(outputPath,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
      status=poc_def_att(outputPath,poc_att_t("history",cmd));
      }
    
    }
  else while(true){
    
    double *elevation=0;
    double *masses;
    static T *deformation=0,*displacement=0,*east=0,*north=0,*pot=0,*lsa=0;
    static float *surfacePressure=0;
    static poc_var_t var;
    
    if(frame<=0)
      unlink(outputPath);
    
    if(zone!=0){
      static poc_data_t<double> input;
      
      if(frame<=0){
        status=input.init(zone,elevAmplitudeName);
        if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_data_t::init(\"%s\",\"%s\") error", zone,elevAmplitudeName);
//         #warning COMMENT LINE BELOW
//         updatemin(&input.nframes,14);
        }
      else if(frame>=input.nframes)
        TRAP_ERR_EXIT(0,"end of computation");
      
      status=input.read_data(zone,frame);
      if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"poc_data_t::read_data(\"%s\",%d) error", zone,frame);
      
      elevation=toP0("elevation",inputMesh,discretisation,input.data);
      
      if(frame<=0 and debugPath!=0 and inputMesh!=0){
        /* using another path as the input mesh may be different from the output mesh */
        string path=replace(debugPath,".nc","-0.nc");
        unlink(path.c_str());
        
        int *P0_ordinals=new int[inputMesh->ntriangles];
        for(t=0;t<inputMesh->ntriangles;t++)
          P0_ordinals[t]=t;
        
        status=poc_put_mesh_vara(path.c_str(),*inputMesh,LGP0,"P0_ordinals",-1,P0_ordinals,"",0);
        NC_CHKERR_BASE_LINE(status,"poc_put_mesh_vara(\""+path+"\",,,\"P0_ordinals\",...) error");
        
        delete[] P0_ordinals;
        
        status=poc_put_mesh_vara(path.c_str(),*inputMesh,LGP0,"eta_P0",-1,elevation,"m",1);
        NC_CHKERR_BASE_LINE(status,"poc_put_mesh_vara(\""+path+"\",,,\"a_eta_P0\",...) error");
        
        status=poc_def_att(path,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
        status=poc_def_att(path,poc_att_t("history",cmd));
        }
      }
    
    static poc_data_t<double> times;
    double time;
    
    if(outputMesh!=0 and inputMesh==0){
      if(frame>0) TRAP_ERR_EXIT(-1,"UNIFORM 1m INPUT not relevant for time series\n");
      STDOUT_BASE_LINE("FAKING UNIFORM 1m INPUT\n");
      elevation=aset(outputMesh->ntriangles,1.);
      }
    else if(outputMesh!=0 and inputMesh!=0 and outputMesh!=inputMesh){
      STDOUT_BASE_LINE("reinterpolating\n");
      int last=-1;
      
      double *oldElevation=elevation;
      if(uTide!=0 or vTide!=0){
        TRAP_ERR_EXIT(ENOEXEC,"angular momentum not coded yet for timeseries\n");
        STDERR_BASE_LINE("WARNING: SPEED REINTERPOLATION NOT CODED YET: CANCELING ANGULAR MOMENTUM.");
        deletep(&uTide);
        deletep(&vTide);
        deletep(&depth);
        }
      elevation=new double[outputMesh->ntriangles];
      const double mask=0.;/* MUST BE 0. */
      
      static sequence_t< UGfield_t<float> > surfacePressureStream;
      
      if(inputPath!=0 and SP_ampName!=""){
        
        if(frame<=0){
          status=times.init(zone,"time",0,"");
          status=times.read_data(zone);
          }
        
        time=times.data[frame];
        
        if(frame<=0){
          discretisation_t *descriptor;
          status=discretisation_init(outputMesh, LGP0, 0);
          descriptor=get_descriptor_address(*outputMesh, LGP0);
          surfacePressureStream=initialise_UG_sequence(time,"",inputPath,SP_ampName,2,outputMesh,descriptor,filestream2UGfield_processor,0);
          surfacePressure=new float[descriptor->nnodes];
          }
        
        status=surfacePressureStream.interpolate(time,surfacePressure);
        if(status!=NC_NOERR) NC_TRAP_ERROR(return,status,1,"surfacePressureStream.interpolate(%g,) error",time);
        status=poc_put_vara(outputPath,times.info,frame,&time);
        }
      
      for(t=0;t<outputMesh->ntriangles;t++){
        double x,y; /*< interpolation coordinates */
        double interpolated;
        int verbose=0;
        
        fe_get_triangle_lon_lat(*outputMesh,t,&x,&y);
        
        //verbose= t==444196;
        if(verbose)
          STDERR_BASE_LINE("[%d](%g;%g):",t,x,y);
        
        n=fe_detect(*inputMesh,inputMesh->triangles,list,listGrid,x,y,last);
        status=fe_interpolate2D(*inputMesh,LGP0,oldElevation,mask,x,y,n,&interpolated);
        if(verbose) fprintf(stderr,"[%d] %d %g\n",n,status,interpolated);
        
        if(isnan(interpolated)) TRAP_ERR_EXIT(ENOEXEC,"NAN found at [%d](%g,%g)\n",t,x,y);
        
        if(surfacePressure!=0){
          float *z=&surfacePressure[t];
          if(verbose) fprintf(stderr,"[%d] %g\n",t,status,*z);
          if(*z==surfacePressureStream.frames[0].mask or isnan(*z))
            TRAP_ERR_EXIT(ENOEXEC,"%g found at [%d](%g,%g)\n",z,t,x,y);
          (*z)*=1./(1025*9.81);
          interpolated+=*z;
          }
        
        elevation[t]=interpolated;
        }
      
      delete[]oldElevation;
      }
    else TRAP_ERR_EXIT(ENOEXEC,"programming error: should have exited before\n");
    
    if(surfacePressure!=0)
      status=poc_put_mesh_vara(outputPath,*outputMesh,LGP0,"Pa",frame,surfacePressure,"m",1);
    var.init("P",NC_FLOAT,"sea surface height minus inverse barometer","m");
    status=poc_put_mesh_vara(outputPath,*outputMesh,LGP1,&var,frame,elevation,1);
    
    /* only here because poc_put_mesh_vara() calls nc_create() */
    status=poc_put_vara(outputPath,times.info,frame,&time);
    
    /** \bug 2012-05-22 Damien Allain : using a constant volumic mass */
    const double rho_water=1025.;/*< volumic mass of sea water */
    double waterColumMass;
    
    if(uTide!=0 and vTide!=0 and depth!=0){
      for(t=0;t<outputMesh->ntriangles;t++){
      /* local momentum */
        waterColumMass=area[t]*rho_water*depth[t];
        uTide[t]*=waterColumMass;
        vTide[t]*=waterColumMass;
        }
      }
    
    if(frame<=0 and debugPath!=0){
/*------------------------------------------------------------------------------
      save output mesh */
      unlink(debugPath);
      
      int *P0_ordinals=new int[outputMesh->ntriangles];
      for(t=0;t<outputMesh->ntriangles;t++)
        P0_ordinals[t]=t;
      
      status=poc_put_mesh_vara(debugPath,*outputMesh,LGP0,"P0_ordinals",-1,P0_ordinals,"",0);
      NC_CHKERR_BASE_LINE(status,"poc_put_mesh_vara(\"%s\",,,\"P0_ordinals\",...) error",debugPath);
      
      delete[] P0_ordinals;
      
      status=poc_put_mesh_vara(debugPath,*outputMesh,LGP0,"area",-1,area,"m2",1);
      status=poc_put_mesh_vara(debugPath,*outputMesh,LGP0,"perimeter_area",-1,perimeter_area,"m",1);
      status=poc_def_att(debugPath,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
      status=poc_def_att(debugPath,poc_att_t("history",cmd));
      }
    
/*------------------------------------------------------------------------------
    complex masses */
    masses=massFromAreaAndElevation<double>(outputMesh->ntriangles,area,elevation);
    delete[]elevation;
    
/*------------------------------------------------------------------------------
    call loading() */
    
    if(frame<=0){
      if(not lsaOnly){
        displacement=new T[outputMesh->nvtxs];
        //deformation=new T[outputMesh->nvtxs];
#if HORIZONTAL_COLUMN >0
        east=new T[outputMesh->nvtxs];
        north=new T[outputMesh->nvtxs];
#endif
        pot=new T[outputMesh->nvtxs];
        }
      lsa=new T[outputMesh->nvtxs];
      }
    
    if(lsa!=0)
      loading<vector3_t,matrix3x3_t>(outputMesh->ntriangles,lat,lon,masses,outputMesh->nvtxs,outputLat,outputLon,deformation,displacement,east,north,pot,lsa,outputMesh);
    
/*------------------------------------------------------------------------------
    save results */
    
    if(not lsaOnly){
      var.init("UP",NC_FLOAT,"upward displacement","m");
      status=poc_put_mesh_vara(outputPath,*outputMesh,LGP1,&var,frame,displacement,1);
      if(deformation!=0)
        status=poc_put_mesh_vara(outputPath,*outputMesh,LGP1,"UP_CF",frame,deformation,"m",1);
#if HORIZONTAL_COLUMN >0
      status=poc_put_mesh_vara(outputPath,*outputMesh,LGP1,"EAST",frame,east,"m",1);
      status=poc_put_mesh_vara(outputPath,*outputMesh,LGP1,"NORTH",frame,north,"m",1);
#endif
      var.init("potential",NC_FLOAT,"gravitational potential","m");
      status=poc_put_mesh_vara(outputPath,*outputMesh,LGP1,&var,frame,pot,1);
      }
    if(lsa!=0){
      var.init("LSA",NC_FLOAT,"loading_and_self-attraction","m");
      status=poc_put_mesh_vara(outputPath,*outputMesh,LGP1,&var,frame,lsa,1);
      }
    
    if(frame<=0){
      status=poc_def_att(outputPath,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
      status=poc_def_att(outputPath,poc_att_t("history",cmd));
      }
    
    delete[]masses;
    frame++;
    }
  
  outputMesh->destroy();
  
  deletep(&area);
  
  return status;
}
#undef T
#undef NC_T
#undef NC_FILL_T
#undef TS


#define T double
#define NC_T NC_DOUBLE
#define NC_FILL_T NC_FILL_DOUBLE
#define TS "double"
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int loading_RG2UG()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// tidal loading with an unstructured grid (regular grid input)
/*----------------------------------------------------------------------------*/
{
STDERR_BASE_LINE("using " TS "\n");
  int t,c,n;//indexes : triangles, corners, nodes
  vertex_t *node;
  triangle_t *triangle;
  double x,y;//interpolation coordinates
  double *lon,*lat;//triangle coordinates
  double *outputLon,*outputLat;//node coordinates

  int status;
  grid_t inputGrid;
  mesh_t mesh;
  double *area;
  complex<float> *inputElevation,inputMask,interpolated;
  complex<double> *elevation;
  int nMasked;
  complex<double> *masses;
  complex<T> *deformation=0,*displacement,*east=0,*north=0,*pot,*lsa,mask;

  if(topoPath!=NULL)STDERR_BASE_LINE("ignoring topography from %s\n",topoPath);
  
  /// <h1>It loads the elevations</h1> with load_atlas()
  status=load_atlas(inputPath, elevAmplitudeName, elevPhaseName, &inputGrid,&inputElevation,&inputMask,POLAR,0);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"load_atlas error with %s and %s in %s", elevAmplitudeName, elevPhaseName,inputPath);
  
  /// <h1>It loads the unstructured grid</h1> with fe_readgeometry()
  status=fe_readgeometry(zone,&mesh);
  if(status!=NC_NOERR)NC_TRAP_ERROR(return,status,1,"fe_readgeometry error with %s",zone);
  /// <h2>It interpolates the elevations to the center of the triangles</h2>
  elevation=new complex<double>[mesh.ntriangles];
  area=new double[mesh.ntriangles];
  lat=new double[mesh.ntriangles];
  lon=new double[mesh.ntriangles];
  nMasked=0;
  for(t=0;t<mesh.ntriangles;t++){
    triangle=&mesh.triangles[t];
    y=0.;
    x=0.;
    n=triangle->vertex[0];
    node=&mesh.vertices[n];
    for(c=0;c<3;c++){
      n=triangle->vertex[c];
      y+=geo_recale(mesh.vertices[n].lat,node->lat);
      x+=geo_recale(mesh.vertices[n].lon,node->lon);
      }
    y/=3.;
    x=map_recale(inputGrid,x/3);
    lat[t]=y;
    lon[t]=x;
    status=map_interpolation(inputGrid, inputElevation, inputMask, x,y,&interpolated);
    if(abs(x<1) && abs(y<1)){
      //STDERR_BASE_LINE("%gE %gN %gm\n",x,y,abs(interpolated));
      }
    if(interpolated==inputMask){
      nMasked++;
      area[t]=0.;
      elevation[t]=0.;
      continue;
      }
    elevation[t]=(complex<double>)interpolated;
    area[t]=triangle->TrueArea;
    }
  STDERR_BASE_LINE("%d/%d masked\n",nMasked,mesh.ntriangles);
  delete[]inputElevation;
  if(debugPath!=0){
    unlink(debugPath);
    status=archiving_UGdummy2D(debugPath,mesh,"Ha","Hg","m",elevation,mask,LGP0);
    }
  
/*---------------------------------------------------------------------*//**<h1>
  calculates complex masses with mass() </h1>*/
  masses=massFromAreaAndElevation<complex<T>  >(mesh.ntriangles,area,elevation);
  delete[]area;
  delete[]elevation;

/*---------------------------------------------------------------------*//**<h1>
  finally calls loading() </h1>*/
  outputLat=new double[mesh.nvtxs];
  outputLon=new double[mesh.nvtxs];
  for(n=0;n<mesh.nvtxs;n++){
    node=&mesh.vertices[n];
    outputLat[n]=node->lat;
    outputLon[n]=node->lon;
    }
//   deformation=new complex<T>[mesh.nvtxs];
  displacement=new complex<T>[mesh.nvtxs];
#if HORIZONTAL_COLUMN >0
  east=new complex<T>[mesh.nvtxs];
  north=new complex<T>[mesh.nvtxs];
#endif
  pot=new complex<T>[mesh.nvtxs];
  lsa=new complex<T>[mesh.nvtxs];
  
  unlink(outputPath);
  
  loading<cvector3_t,zmatrix3x3_t>(mesh.ntriangles,lat,lon,masses,mesh.nvtxs,outputLat,outputLon,deformation,displacement,east,north,pot,lsa,&mesh);
  
/*---------------------------------------------------------------------*//**<h1>
  Then it saves the results with archiving_UGdummy2D() </h1>*/
  status=archiving_UGdummy2D(outputPath,mesh,"UP_a","UP_G","m",displacement,mask,LGP1);
  if(deformation!=0)
    status=archiving_UGdummy2D(outputPath,mesh,"UP_CF_a","UP_CF_G","m",deformation,mask,LGP1);
#if HORIZONTAL_COLUMN >0
  status=archiving_UGdummy2D(outputPath,mesh,"EAST_a","EAST_G","m",east,mask,LGP1);
  status=archiving_UGdummy2D(outputPath,mesh,"NORTH_a","NORTH_G","m",north,mask,LGP1);
#endif
  status=archiving_UGdummy2D(outputPath,mesh,"potential_a","potential_G","m",pot,mask,LGP1);
  status=archiving_UGdummy2D(outputPath,mesh,"LSA_a","LSA_G","m",lsa,mask,LGP1);
  status=poc_def_att(outputPath,poc_att_t("production","constructed around " __LINE_FILE_PACKAGE_REVISION));
  status=poc_def_att(outputPath,poc_att_t("history",cmd));
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void loading_RG()

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// tidal loading with a regular grid
/*----------------------------------------------------------------------------*/
{
STDERR_BASE_LINE("using " TS "\n");
  int n,i,j;//indexes : node, columns, rows
  double x,y;//interpolation coordinates
  bool debug=false;

  int status;
  grid_t inputGrid,grid,outputGrid;
  int gridn,outputGridn;
  double *area;
  complex<T> *inputElevation=0,elevationMask,*elevation;
  complex<T> *masses;
  complex<double> total;
  complex<T> *deformation=0,*displacement,*east=0,*north=0,*pot,*lsa;
  
/*---------------------------------------------------------------------*//**<h1>
  loads the elevations </h1>*/
  if(inputPath[0]){
    poc_cdata_t inputData;
    status=inputData.init(inputPath, elevAmplitudeName, elevPhaseName);
    status=poc_get_grid(inputPath, inputData.info, &inputGrid);
    elevationMask=inputData.mask;
    swapValues(&inputElevation,&inputData.data);
    }
  else{///or fakes a test set
    //x=1./32;map_set2Dgrid(&inputGrid,-3.+x/2,49.+x/2,-2.,50.,x);
    x=2.;map_set2Dgrid(&inputGrid,x/2,-90.+x/2,360.,90.,x);
    status=map_completegridaxis(&inputGrid,2);
    STDERR_BASE_LINE("FAKE inputGrid:");inputGrid.print();
    free(zone);zone=NULL;
    free(topoPath);topoPath=NULL;
    gridn=inputGrid.nx*inputGrid.ny;
    inputElevation=new complex<T>[gridn];
    for(n=0;n<gridn;n++){
      inputElevation[n]=1.;
      }
    }
  if(zone!=NULL){/// <h2>If the zone is specified</h2>
    ///it sets the regular grid accordingly
    grid=get_zonegrid_shifted(zone);
    status=map_completegridaxis(&grid,2);
    STDERR_BASE_LINE("grid:");grid.print();
    ///and it interpolates the elevations to the required resolution
    elevation=new complex<T>[grid.nx*grid.ny];
    for (j=0;j<grid.ny;j++) {
      for (i=0;i<grid.nx;i++) {
        n=i+grid.nx*j;
        y=grid.y[n];
        x=map_recale(inputGrid,grid.x[n]);
        status=map_interpolation(inputGrid, inputElevation, elevationMask, x,y,&elevation[n]);
        }
      }
    delete[]inputElevation;
    
    if(debugPath){
      STDERR_BASE_LINE("saving interpolated elevation\n");
      poc_var_t var;
      var=poc_save_grid(debugPath, grid, __FILE__, __LINE__);
      var.init("elevation",NC_DOUBLE,"sea_surface_height_due_to_non_equilibrium_ocean_tide","",NC_FILL_DOUBLE);
      status=poc_put_cvara(debugPath, var, 0, elevation);
      STDERR_BASE_LINE("status=%d\n",status);
      }
    }
  else{/// <h2>Otherwise</h2>
    ///it takes the same zone as the input.
    grid=inputGrid;
    STDERR_BASE_LINE("(input)grid:");grid.print();
    elevation=inputElevation;
    }
  gridn=grid.Hsize();
  if(topoPath!=NULL){/// <h2>If the bathymetry is specified</h2>
    short int *topo,topoMask;
    float topoValue;
    grid_t topoGrid;
    status=topo_loadfield(topoPath, &topoGrid, &topo, &topoMask, debug);
    ///it masks the elevations accordingly
    for (j=0;j<grid.ny;j++) {
      for (i=0;i<grid.nx;i++) {
        n=i+grid.nx*j;
        y=grid.y[n];
        x=map_recale(topoGrid,grid.x[n]);
        status=map_interpolation(topoGrid, topo, topoMask, x,y,&topoValue);
        if(topoValue>=0){
          elevation[n]=elevationMask;
          }
        }
      }
    delete[]topo;
    STDERR_BASE_LINE("saving masked elevation\n");
    poc_var_t var;
    var=poc_save_grid("maskedLoadingElevation.nc", grid, __FILE__, __LINE__);
    var.init("elevation",NC_DOUBLE,"sea_surface_height_due_to_non_equilibrium_ocean_tide","",NC_FILL_DOUBLE);
    status=poc_put_cvara("maskedLoadingElevation.nc", var, 0, elevation);
    STDERR_BASE_LINE("status=%d\n",status);
    }
  
/*---------------------------------------------------------------------*//**<h1>
  calculates areas and complex masses with mass() </h1>*/
  status=map_completegridaxis(&grid,2);
  
  area=new double[gridn];
  i=0;
  
  for(n=0;n<gridn;n++){
    
    /* skip longitude repeats */
    if(grid.modeH>=0){
      i++;
      if(i==grid.nx){
        i=0;
        }
      }
    else{
      TRAP_ERR_EXIT(ENOEXEC,"Not coded yet for grid.modeH=%d<0\n",grid.modeH);
      }
    
    if(elevation[n]!=elevationMask && i!=0){
      /* integrating */
      //area[n]=(sin((grid.y[n]+grid.dy/2)*d2r)-sin((grid.y[n]-grid.dy/2)*d2r))*grid.dx*d2r*square(MeanEarthRadius);
      /* is not better than this */
      area[n]=cos(grid.y[n]*d2r)*grid.dx*grid.dy*square(d2m);
      continue;
      }
    
    area[n]=0.;
    elevation[n]=0.;
    }
  
  masses=massFromAreaAndElevation<complex<T> >(gridn,area,elevation);
  delete[]area;
  delete[]elevation;
  
/*---------------------------------------------------------------------*//**<h1>
  finally calls loading() </h1>*/
  const int nt=742114,nn=389737;
  outputGrid=grid;
  STDERR_BASE_LINE("outputGrid:");outputGrid.print();
  outputGridn=outputGrid.Hsize();
  
//   deformation=new complex<T>[outputGridn];
  displacement=new complex<T>[outputGridn];
#if HORIZONTAL_COLUMN >0
  east=new complex<T>[outputGridn];
  north=new complex<T>[outputGridn];
#endif
  pot=new complex<T>[outputGridn];
  lsa=new complex<T>[outputGridn];
  
  double cct;
  cct=loading<cvector3_t,zmatrix3x3_t>(gridn,grid.y,grid.x,masses,outputGridn,outputGrid.y,outputGrid.x,deformation,displacement,east,north,pot,lsa,0,&inputGrid);
  cct*=(double)nt/gridn;
  STDOUT_BASE_LINE("This means %gs for the same output with %d inputs and %gs for a %d-point output\n",cct,nt,cct*(double)nn/outputGridn,nn);
  
  /*
  ///It checks volume conservation
  total=0.;
  for(n=0;n<gridn;n++){
    total+=(complex<double>)load[n]*cos(outputGrid.y[n]*d2r)*square(d2m);
    }
  STDOUT_BASE_LINE("Earth compression : %gm3=>%gm\n",abs(total),abs(total)/(4*M_PI*square(MeanEarthRadius)*.63));
  */
  
/*------------------------------------------------------------------------------
  save outputs */
  STDOUT_BASE_LINE("saving to %s: ",outputPath);fflush(stdout);
  poc_var_t var;
  status=poc_save_grid(outputPath,&var,outputGrid,1,0);
  
  var.init("UP",NC_FLOAT);
  printf(var.name+", ");fflush(stdout);
  status=poc_put_cvara(outputPath,var,0,displacement);
  if(deformation!=0){
    var.name="deformation";
    printf(var.name+", ");fflush(stdout);
    status=poc_put_cvara(outputPath,var,0,deformation);
    }
#if HORIZONTAL_COLUMN > 0
  var.name="EAST";
  printf(var.name+", ");fflush(stdout);
  status=poc_put_cvara(outputPath,var,0,east);
  var.name="NORTH";
  printf(var.name+", ");fflush(stdout);
  status=poc_put_cvara(outputPath,var,0,north);
#endif
  var.name="potential";
  printf(var.name+", ");fflush(stdout);
  status=poc_put_cvara(outputPath,var,0,pot);
  var.name="LSA";
  printf(var.name+", ");fflush(stdout);
  status=poc_put_cvara(outputPath,var,0,lsa);
  
  printf("done.\n");fflush(stdout);
}
#undef T
#undef NC_T
#undef NC_FILL_T
#undef TS


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// function that will be ran when the executable is started
/** See the print_help <a href=#func-members>function</a> for its use.
*/
/*----------------------------------------------------------------------------*/
{
  int n;//argument index
  char *keyword;
  
  struct timeval mainbefore;
  cmd=fct_echo( argc, argv);
  STDERR_BASE_LINE("PID: %d\n",getpid());
  gettimeofday(&mainbefore);
  
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if( strcmp(keyword,"-h")==0 or strcmp(keyword,"--help")==0 ) {
          print_help(argv[0]);
          return 0;
          }
        else if( strcmp(keyword,"--lsa-only")==0 ){
#if POT_UP_TO_LSA == 1
          TRAP_ERR_EXIT(ENOEXEC,"not coded for --lsa-only."
            "Please see the definitionS of POT_UP_TO_LSA in the code\n");
#endif
          lsaOnly=true;
          n++;
          }
        else switch (keyword[1]) {
        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshPath= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'r' :
          inputPath= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'v' :
          fprintf(stderr,"*** Use of -v is osolete. Use -e instead ***\n");
          print_help(argv[0]);
          exit(-1);

        case 'e' :
          elevAmplitudeName= argv[n+1];
          elevPhaseName= argv[n+2];
          n++;
          n+=2;
          break;

        case 'u' :
          uAmpName= argv[n+1];
          uPhaName= argv[n+2];
          vAmpName= argv[n+3];
          vPhaName= argv[n+4];
          n++;
          n+=4;
          break;

        case 'o' :
          outputPath= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'D' :
          debugPath= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          topoPath= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          STDOUT_BASE_LINE("unknown option %s\n",keyword);
          print_help(argv[0]);
          exit(-1);
        }
        break;

      default:
        STDOUT_BASE_LINE("unknown option %s\n",keyword);
        print_help(argv[0]);
        exit(-1);
        break;
      }
    free(keyword);
    }
  
  do{
    const bool
      zoneIsFile= zone!=NULL and access(zone,R_OK)==0,
      zoneIsZone= zone!=NULL and not zoneIsFile;
      
    
    if(zoneIsFile or meshPath!=0){
      if(zoneIsFile){
        STDERR_BASE_LINE("The zone is a readable file\n");
        copyTo_tmpdir(&zone,1);
        }
      if(meshPath)
        STDERR_BASE_LINE("continental mesh is: %s\n",meshPath);
      
      if(meshPath==0 and inputPath!=NULL and access(inputPath,R_OK)==0){
        /* If the regular input is readable */
        loading_RG2UG();
        }
      else{
        /* If the regular input is NOT readable */
        loading_UG();
        }
      break;
      }
    
    if(zoneIsZone)STDERR_BASE_LINE("The zone is a name or a resolution : %s\n",zone);
    else STDERR_BASE_LINE("The zone is not specified\n");
    
    if(inputPath==NULL){
      STDOUT_BASE_LINE("*** Missing elevation or zone ***\n");
      print_help(argv[0]);
      exit(-1);
      }
    
/*---------------------------------------------------------------------*//**<h1>
    If the zone is missing or is a grid name or a resolution </h1>*/
    ///calls loading_RG()
    loading_RG();
    
    }while(0);
  
  STDOUT_BASE_LINE(" ^^^^^^^^^^^^^ end of computation, that took %#g s ^^^^^^^^^^^^^\n",difftime(mainbefore));
  exit(0);
}
