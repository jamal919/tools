
/**************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

***************************************************************************/
/** \file

\brief tidal operation prototypes

Declares functions defined in tides.cpp tides01.cpp tides02.cpp harmonic-01.cpp harmonic-02.cpp
\sa tides.def
*/
/*----------------------------------------------------------------------------*/

#ifndef _TIDES_H
#define _TIDES_H

#include "constants.h"
#include "tools-define.h"

#include "poc-netcdf-data.hpp"

const double pi2 = M_PI*2.;

const double dph2rps=d2r/3600.;
const double dph2rpd=24.*d2r;
const double dph2cpd=24./360.;

/** See http://hpiers.obspm.fr/eop-pc/models/constants.html */

const double sy2d=365.256363; ///< number of days per sidereal year

const double d2s=86400.; ///< number of seconds per day

const double jc2d=36525.;///< number of days per Julian century

/* NOTE: the Earth spins sy2d+1 times per year ! */
const double two_Omega = 2*2*M_PI/(d2s*sy2d/(sy2d+1));///<2*Omega, in rad/s
  //two_Omega = 1.45444104e-4;//FALSE value = 2*2*M_PI/d2s

#include "tools-structures.h"
#include "netcdf-classes.h"
#include "poc-time.h"
#include "functions.h"
#include "harmonic.h"

#include "ascii.h"

#define POLAR     0
#define CARTESIAN 1

#define IHO_LIST_URL "http://www.iho.int/mtg_docs/com_wg/IHOTC/IHOTC_Misc/TWLWG_Constituent_list.pdf"
#define WAVELISTURL "http://www.archive.org/details/manualofharmonic00schu"
#define WAVELISTPARTIALREF "Schureman (1940, Tab2 pp.164-167)"
#define WAVELISTREF WAVELISTPARTIALREF " " WAVELISTURL

#define MAXIMUM_WAVE_NAME_LENGTH 10

class compound_tools_t {
private :
public :
  tidal_wave w[10];             /*  specify generating waves                           */
  int c[10];                    /*                                                     */
  int n;                        /*                                                     */
  char name[20];                /*                                                     */
  
  compound_tools_t() {
    n=0;
    }

  compound_tools_t init(const char *wname) {
    strcpy(name,wname);
    n=0;
    return *this;
    }
  ///\todo 2011-10-10 Damien Allain : automatic calls of add() deduced from the name of the wave ?
  compound_tools_t add(tidal_wave wave, int coefficient) {
    w[n]=wave;
    c[n]=coefficient;
    n++;
    return *this;
    }

  tidal_wave realize(){
    int k;
    tidal_wave wave;
    wave.reset();
    for(k=0;k<n;k++) {
      wave.nT+=c[k]*w[k].nT;
      wave.ns+=c[k]*w[k].ns;
      wave.nh+=c[k]*w[k].nh;
      wave.np+=c[k]*w[k].np;
      wave.np1+=c[k]*w[k].np1;
      wave.nksi+=c[k]*w[k].nksi;
      wave.nnu0+=c[k]*w[k].nnu0;
      wave.nnu1+=c[k]*w[k].nnu1;
      wave.nnu2+=c[k]*w[k].nnu2;
      wave.Qu+=c[k]*w[k].Qu;
      wave.Ra+=c[k]*w[k].Ra;
      wave.shift+=c[k]*w[k].shift;
      wave.code=0;
      }
    strcpy(wave.name,name);
    wave.init("s");
    return(wave);
    }

};

class compound_t {
private :
public :
/*-----------------------------------------------------------------------------------------------------
  physical parameters*/
  tidal_wave  wave;
  //tide2D_t    *state;
  spectrum_t  generating;
  vector<string> formula;
  int id,n;

  compound_t() {
    n=0;
//    formula=string("");
    }
};

class astro_angles_t{
public:
  double t_julian;///<time elapsed between 1900-01-01 00:00 and the reference in julian centuries (of 36525 days)
  ///\note 1900-01-01 00:00  is the time reference for astronomical abaques
  /* all that follow are calculated in astronomic_angle() */
  double sh_T,sh_h,sh_s,sh_p1,sh_p;
  double sh_xi,sh_nu,sh_R,sh_x1ra,sh_nuprim,sh_nusec,sh_Qu;
  double sh_I,sh_N;
};


class tidaldata_t {
private :
public :
  double lon, lat;       /*  location of the data   */
  float  *Ha, *Hg;       /*  specified tidal data   */
  float  *Ua, *Ug;       /*  specified tidal data   */
  float  *Va, *Vg;       /*  specified tidal data   */
  
  tidaldata_t() {
    Ha=Hg=Ua=Ug=Va=Vg=0;
    }
  void destroy() {
    deletep(&Ha);
    deletep(&Hg);
    deletep(&Ua);
    deletep(&Ug);
    deletep(&Va);
    deletep(&Vg);
    }
} ;


class parameter_t {
private :
public :
/*----------------------------------------------------------------------------------------------------
  geometry parameters*/
  double *h, *dhdx, *dhdy;          /* mean depth, depth gradients                                  */
  double *h__;                      /* mean depth, depth gradients                                  */
  double **sigma;                   /* sigma coordinates                                            */
  double *lon, *lat;                /* longitude, latitude                                          */
  double *C, *S;                    /* cosine and sine of the latitude                              */
  float  *lmm;                      /* lumped mass matrix                                           */
/*----------------------------------------------------------------------------------------------------
  geometry parameters*/
  int *countdown, *unstable;        /*                                                               */
  int *sponge;                      /*                                                               */
  int *uCode,   *zCode;             /*                                                               */
  int *uParent, *zParent;           /*                                                               */
  int *zAncestor;                   /* subcycle transversal information                              */
/*-----------------------------------------------------------------------------------------------------
  physical parameters*/
  float  *vwr,*z0, *z03D;           /* wave drag rugosty, bottom roughness 2D and 3D                 */
  float  *rlinear, **rlinear3D;     /* linear friction coeffcient, 2D, 3D                            */
  float  *Cd,*Cd3D;                 /* quadratic friction coefficient, 2D, 3D                        */
  float  *FrictionRatio;            /*                                                               */
  double **Kv3D,**epsilon,**tke;    /* vertical diffusivity, dissipation rate, turbulent kinetic energy */
  float  *wdc1, *Nbar, *celerity;   /* wave drag coeffcient (slope effect)                           */
  double *u0, *u0__;                /* background velocity (for friction)                            */
/*-----------------------------------------------------------------------------------------------------
  tidal parameters*/
  hconstant_t *Htide,*Utide,*Vtide; /*  Computed tidal constants                                     */
  hconstant_t *LSA;                 /*  Loading and self-attraction                                  */
  int    *have_LSA;                 /*                                                               */
  int    nLSA;                      /*                                                               */
  int    colocation;                /* flag indicating pressure and velocity colocation              */

  parameter_t() {
    h=0; dhdx=0; dhdy=0;
    h__=0;
    sigma=0;
    lon=0; lat=0;
    C=0; S=0;
    lmm=0;
    countdown=0; unstable=0;
    sponge=0;
    uCode=0; zCode=0;
    uParent=0; zParent=0;
    zAncestor=0;
    vwr=0; z0=0;
    rlinear=0; rlinear3D=0; Cd=0; Cd3D=0;
    FrictionRatio=0;
    wdc1=Nbar=celerity=0;
    u0=0; u0__=0;
    Htide=0; Utide=0; Vtide=0;
    have_LSA=0;
    LSA=0;
    nLSA=-1;
    colocation=-1;
    Cd3D=0;
    Kv3D=0;
    epsilon=0;
    tke=0;
    }

  void import_z(int target, int source, parameter_t data) {
    lon[target]=data.lon[source];
    lat[target]=data.lat[source];
    h[target]=data.h[source];
//     dhdx[target]=data.dhdx[source];
//     dhdy[target]=data.dhdy[source];
    if(data.sigma !=0) sigma[target]=data.sigma[source];
    C[target]=data.C[source];
    S[target]=data.S[source];
    lmm[target]=data.lmm[source];
    zCode[target]=data.zCode[source];
    zParent[target]=data.zParent[source];
    zAncestor[target]=data.zAncestor[source];
//     countdown[target]=data.countdown[source];
//     unstable[target]=data.unstable[source];
    sponge[target]=data.sponge[source];
    u0__[target]=data.u0__[source];
    Htide[target]=data.Htide[source];
    Utide[target]=data.Utide[source];
    Vtide[target]=data.Vtide[source];
    if(data.Cd3D!=0) Cd3D[target]=data.Cd3D[source];
    //Kv3D[target]=data.Kv3D[source];
    if(data.LSA!=0) LSA[target]=data.LSA[source];
    have_LSA=data.have_LSA;
    switch(colocation) {
      case 0:
        u0__[target]=data.u0__[source];
        break;

      case 1:
        break;

      default:
        printf("colocation flag %d not valid\n",colocation);
        exit(-1);
        break;
      }
    }

  void import_u(int target, int source, parameter_t data) {
    h__[target]=data.h__[source];
    dhdx[target]=data.dhdx[source];
    dhdy[target]=data.dhdy[source];
    uParent[target]=data.uParent[source];
    uCode[target]=data.uCode[source];
    vwr[target]=data.vwr[source];
    rlinear[target]=data.rlinear[source];
    Cd[target]=data.Cd[source];
    z0[target]=data.z0[source];
    FrictionRatio[target]=data.FrictionRatio[source];
    wdc1[target]=data.wdc1[source];
    Nbar[target]=data.Nbar[source];
    celerity[target]=data.celerity[source];
    u0[target]=data.u0[source];
    }


  void destroy() {
    int n;
    deletep(&h);
    deletep(&dhdx);
    deletep(&dhdy);
    deletep(&lon);
    deletep(&lat);
    deletep(&C);
    deletep(&S);
    deletep(&lmm);
    deletep(&countdown);
    deletep(&unstable);
    deletep(&sponge);
    deletep(&zCode);
    deletep(&zParent);
    deletep(&zAncestor);
    deletep(&uCode);
    deletep(&uParent);
    deletep(&vwr);
    deletep(&rlinear);
    deletep(&Cd);
    deletep(&u0);
    deletep(&z0);
    deletep(&FrictionRatio);
    deletep(&wdc1);deletep(&Nbar);deletep(&celerity);
    deletep(&Htide);
    deletep(&Utide);
    deletep(&Vtide);
    deletep(&Cd3D);
    deletep2D(&Kv3D,-1);
    deletep2D(&epsilon,-1);
    deletep2D(&tke,-1);
    deletep2D(&sigma,-1);
    if(LSA!=0) {
      for(n=0;n<nLSA;n++) LSA[n].destroy();
      nLSA=0;
      }
    deletep(&LSA);
    deletep(&have_LSA);
    switch(colocation) {
      case 0:
        deletep(&h__);
        deletep(&u0__);
        break;

      case 1:
        break;

      default:
        printf("colocation flag %d not valid\n",colocation);
        exit(-1);
        break;
      }
    }
//   ~parameter_t() {
//     this->destroy();
//     }
};

class tidalOBC_t {
private :
public :
  double lon, lat;                  /* longitude, latitude                   */
/*-----------------------------------------------------------------------------
  tidal parameters*/
  hconstant_t Htide,Utide,Vtide;    /* Computed tidal constants              */
  int node;                         /* node index                            */
  
  tidalOBC_t() {
//    Htide=Utide=Vtide=0;
    }
};



class tseries_t {
private :
public :

  double **x;
  double *t,sampling;
  string time_zone;
  char   time_unit;
  double mask;
  double lon,lat;
  date_t first,last,origin;
  int n,nparam;
  mooring_t mooring;

  tseries_t() {
    x=NULL;
    t=NULL;
    n=0;
    lon=lat=0.;
    nparam=0;
    sampling=-1;
    mask=1.e+10;
    }

  tseries_t(double *time, double **obs, double mask, int nobs, int ncol) {
    n=nobs;
    nparam=ncol;
    t=new double[n];
    for(size_t i=0;i<n;i++) t[i]=time[i];
    x=new double*[ncol];
    for (size_t k=0;k<nparam;k++) {
      x[k]=new double[n];
      for(size_t i=0;i<n;i++) x[k][i]=obs[k][i];
      }
    lon=lat=0.;
    this->mask=mask;
    }

  tseries_t(int ncol) {
    n=0;
    nparam=ncol;
    x=new double*[ncol];
    }

  tseries_t(int nobs, int ncol) {
    n=nobs;
    nparam=ncol;
    t=new double[n];
    x=new double*[nparam];
    for (size_t k=0;k<nparam;k++) {
      x[k]=new double[n];
      }
    }

  void allocate(int nobs, int ncol) {
    n=nobs;
    nparam=ncol;
    t=new double[n];
    x=new double*[nparam];
    for (size_t k=0;k<nparam;k++) {
      x[k]=new double[n];
      }
    }

  void redim(int nadd) {
    double **extended=new double*[nparam+nadd];
    for (size_t k=0;k<nparam;k++) {
      extended[k]=x[k];
      }
    for (size_t k=nparam;k<nparam+nadd;k++) {
      extended[k]=new double[n];
      }
    delete[] x;
    x=extended;
    nparam+=nadd;
    }

  size_t size(int param, double value) const {
    size_t nvalid=0;
    for (size_t k=0;k<n;k++) {
      if(x[param][k]==value) nvalid++;
      }
    return(nvalid);
    }

  double length(char unit='s') const {
    double span;
    span=t[n-1]-t[0];
    return(span);
    }

  size_t maxsize() const {
    size_t s;
    s=this->length()/sampling+1;
    return(s);
    }

  int normalize(int target) const {
    statistic_t s;
    s=get_statistics(this->x[target], this->mask, this->n, 0);
    
    for(int k=0;k<this->n;k++) {
      if(this->x[target][k]==this->mask) continue;
      this->x[target][k]=(this->x[target][k]-s.mean)/s.std;
      }
    
    return(0);
    }

  int split(vector<tseries_t> splitted) {
    
    for (size_t k=0;k<this->nparam;k++) {
      tseries_t *tmp=new tseries_t;
      *tmp=tseries_t(t, &(x[k]), mask, n, 1);
      splitted.push_back(*tmp);
      }
    
    return(0);
    }

  int split(tseries_t* & splitted) {
    splitted=new tseries_t[this->nparam];
    for (size_t k=0;k<this->nparam;k++) {
      splitted[k]=tseries_t(t, &(x[k]), mask, n, 1);
      }
    
    return(0);
    }

  tseries_t reduce(date_t start, date_t end) const {
    date_t reference=date_t(1950,1,1,0.);
    int i, first, last, count;
    
    double t1 = julian_day(start)-julian_day(reference)+(double) start.second/d2s;
    double t2 = julian_day(end)  -julian_day(reference)+(double) end.second/d2s;
    
    count=0;
    i=0;
    while(this->t[i] < t1 - 1e-5) {
      if(i==this->n-1) break;
      i++;
      }
    first=i;
    while(this->t[i]<=t2) {
      if(i>=this->n) break;
      i++;
      count++;
      }
    last=i;
    
    tseries_t *reduced= new tseries_t(count,this->nparam);
    
    reduced->lon=this->lon;
    reduced->lat=this->lat;
    reduced->mask=this->mask;
    reduced->mooring=this->mooring;
    
    for(i=first;i<last;i++) {
      reduced->t[i-first]=this->t[i];
      for (size_t k=0;k<this->nparam;k++) {
        reduced->x[k][i-first]=this->x[k][i];
        }
      }
    return(*reduced);
    }

  tseries_t subsample(int increment) {
    int i, count;
    
    count=this->n/increment+1;
    
    tseries_t *reduced= new tseries_t(count,this->nparam);
     
    reduced->lon=this->lon;
    reduced->lat=this->lat;
    reduced->mask=this->mask;
    reduced->mooring=this->mooring;
     
    count=0;
    for(i=0;i<this->n;i+=increment) {
      reduced->t[count]=this->t[i];
      for (size_t k=0;k<this->nparam;k++) {
        reduced->x[k][count]=this->x[k][i];
        }
      count++;
      }

    reduced->n=count;

    return(*reduced);
    }

  int resample_count(double sampling,double t0=NAN) const {
    int count;
    
    if(isnan(t0))
      t0=this->t[0];
   
    count=(this->t[this->n-1]-t0)/sampling+1;
    
    return(count);
    }

  tseries_t resample(double t0, double sampling, double max_gap) {
    int i, count, i0=-1;
    int status;
   
    count=resample_count(sampling,t0);
    
    tseries_t *resampled= new tseries_t(count,this->nparam);
    
    resampled->lon=this->lon;
    resampled->lat=this->lat;
    resampled->mask=this->mask;
    resampled->mooring=this->mooring;
    
    for(i=0;i<resampled->n;i++) {
      resampled->t[i]=t0+(double) i*sampling;
      for (size_t k=0;k<this->nparam;k++) {
        double z;
        status=map_interpolate1D(this->x[k],this->t, this->mask, this->n, resampled->t[i], &z, max_gap, 0, &i0);
        resampled->x[k][i]=z;
        }
      }
    
    return(*resampled);
    }

  tseries_t resample(double sampling, double max_gap) {
    double t0=this->t[0];
    tseries_t resampled=this->resample(t0, sampling, max_gap);
    return(resampled);
    }

  tseries_t resample(double *time, int count, double max_gap, int use_prior=1) {
    int i, i0_=-1, *i0=0;
    int status;
    
    if(use_prior)
      i0=&i0_;
    
    tseries_t *resampled= new tseries_t(count,this->nparam);
    
    resampled->lon=this->lon;
    resampled->lat=this->lat;
    resampled->mask=this->mask;
    resampled->mooring=this->mooring;
    
    for(i=0;i<resampled->n;i++) {
      resampled->t[i]=time[i];
      for (size_t k=0;k<this->nparam;k++) {
        double z;
        status=map_interpolate1D(this->x[k],this->t, this->mask, this->n, resampled->t[i], &z, max_gap, 0, i0);
        resampled->x[k][i]=z;
        }
      }
    
    return(*resampled);
    }
    
  statistic_t stat(int param){
    statistic_t s;
    s=get_statistics(x[param], mask, n, 1);
    }

  void destroy(){
    deletep2D(&x,nparam);
    deletep(&t);
    x=0;
    t=0;
    n=0;
    lon=lat=0.;
    nparam=0;
    sampling=-1;;
    mask=1.e+10;
    }
  
};

//from tide.cpp
extern  int tide_addwave(spectrum_t *list, tidal_wave wave);

extern  void astronomic_angle(double tj, int verbose) __attribute__((safe_deprecated("Add a local astro_angles_t* as 1st argument to your call")));
extern  void astronomic_angle(astro_angles_t *astro_angles, double tj, int verbose);
extern void print_astro_angles(const astro_angles_t &astro_angles);

extern  double nodal_factor(int formula) __attribute__((safe_deprecated("Add a local astro_angles_t as 1st argument to your call")));
extern  double nodal_factor(const astro_angles_t &astro_angles, int formula);

extern  double greenwhich_argument(tidal_wave w) __attribute__((safe_deprecated("Add a local astro_angles_t as 1st argument to your call")));
extern  double greenwhich_argument(const astro_angles_t &astro_angles, tidal_wave w);

extern  double nodal_phase(tidal_wave w) __attribute__((safe_deprecated("Add a local astro_angles_t as 1st argument to your call")));
extern  double nodal_phase(const astro_angles_t &astro_angles, tidal_wave w);

extern int harmonic_coefficients(int, date_t, double *, spectrum_t, harmonic_t *, int);

extern  double potential(tidal_wave w, float x, float y) __attribute__((safe_deprecated("Add a local astro_angles_t as 1st argument to your call")));
extern  double potential(const astro_angles_t &astro_angles, tidal_wave w, float x, float y);
extern  int tidal_equilibrium(tidal_wave w, double x, double y, float *amp, float *pha);

extern  spectrum_t initialize_tide(astro_angles_t *astro_angles, date_t reference);
extern  spectrum_t initialize_tide(bool Z0=true);

extern  void tide_wavearray_initC(spectrum_t *spectrum_list, int nwave);
extern  int index_from_name(const spectrum_t & s, const char *name);
extern  tidal_wave wave_from_name(const char *name);
extern  int initialize_omega(spectrum_t *s);
extern  void init_argument(astro_angles_t *astro_angles, date_t date, int verbose=0);
extern void init_argument(astro_angles_t *astro_angles, double startd, int verbose=0);
extern  void tides_compound();

extern double tide_alias(tidal_wave w, double Trepet);
extern void spectrum_define_alias(spectrum_t *s, double Trepet, FILE *out=0);

extern void spectrum_redef(spectrum_t *s);

extern void spectrum_isSolvable(spectrum_t s, double duration, int *keep, int *resolve, FILE *out=0);
extern void spectrum_isRedundant(spectrum_t s, int *keep);
extern void spectrum_checkRayleigh(spectrum_t s, double duration, int *keep, int *resolve, FILE *log =0);
extern void spectrum_isAdmittanceSolvable(spectrum_t s, double duration, int *keep, int *resolve, int *deduce, FILE *out=0);

extern void spectrum_reset(spectrum_t& s, const string spectrumType, const double Trepet, FILE *out=0);


extern mgr_data_t *harmonic_analysis_with_parsing(const double *serie, double mask,const double *time, double *residuals, double *mean, int n, spectrum_t prior, spectrum_t & s, int nodal, FILE *log=0);
extern mgr_data_t *harmonic_analysis_without_parsing(double *, double *, double *,  double *, int, spectrum_t, int, int *, int *, FILE *);

extern mgr_data_t *harmonic_analysis(const double *serie, double *residuals,const double *time, int n, const spectrum_t & s, double averaging_time, int* & keep, int* & deduce, const hconstant_t & prior, int nodal, double maxratio, FILE *log, int tag=-1);
extern mgr_data_t *harmonic_analysis(const double *serie, double *residuals,const double *time, int n,const spectrum_t & prior,double averaging_time, double repetition, spectrum_t & s, int nodal, double maxratio, FILE *log);

extern void matrix_member(const double *A, const spectrum_t & s, int k, int l, double *rere, double *imim, double *reim=0, double *imre=0);
extern double matrix_correlation(const double *A, const spectrum_t & s, int k, int l, double rk, double ik);
extern void matrix_auto_correlation(const double *A, const spectrum_t & s, int k, double *re, double *im);
extern void save_full_matrix(const string & path, int *index, int tag, const spectrum_t & s, const double *A, const int *keep=0, const int *deduce=0, double *rhs=0, double *zr=0, double *zi=0, double *error=0);

extern mgr_data_t *check_parsing(double *time, int n, const spectrum_t & s, int nodal, FILE *log);


extern int harmonic_analysis_old(double *,double *,double *, double *,int, spectrum_t);
extern int harmonic_analysis_old(float  *,double *,float  *, float  *,int, spectrum_t);

///Whether the hconstant_t** returned by harmonic_analysis_core() is polar-initialised.
/** If not, it is complex-initialised */
#define harmonic_analysis_constants_USE_polar 0
extern hconstant_t **harmonic_analysis_core(harmonic_t, double, int);

extern int harmonic_predictionTS(double *serie,const double *time, int n,const spectrum_t & s, float *a, float *G, int nodal=1);

extern int harmonic_predictionTS(float  *serie, double *time, int n, spectrum_t s, complex<float> *constants, int nodal=1);
extern int harmonic_predictionTS(double *serie, double *time, int n, spectrum_t s, complex<float> *constants, int nodal=1);
extern int harmonic_predictionTS(double *serie, double *time, int n, mgr_data_t *data, int ndata, int nodal);
extern int harmonic_predictionTS(double *serie, double *time, int n, hconstant_t & data, int nodal);
extern int harmonic_predictionTS(tseries_t &, date_t, date_t, double, hconstant_t &, int);

extern int harmonic_prediction(double *buffer,double time,int n, spectrum_t s, complex<float> **constants, int nodal=1);

extern int harmonic_prediction(const astro_angles_t &astro_angles,double *buffer,double time,int n,const spectrum_t &s,hconstant_t *constants,int nodal=1);

extern void harmonic_init(double *h, double mask, int n, date_t start, double dt, spectrum_t s, harmonic_t *ptr);
extern void harmonic_init_nr(int n, date_t start, double *time, spectrum_t s, harmonic_t *ptr);
//from tides-energy.cpp
extern int harmonic_init01(int nframe, double *time, spectrum_t s, harmonic_t *x, int nodal,const astro_angles_t &astro_angles);
extern int harmonic_init01(int nframe, date_t start, double *time, spectrum_t s, harmonic_t *x, int nodal,const astro_angles_t &astro_angles) __attribute__((safe_deprecated("second argument has no effect")));
//EO from tides-energy.cpp


extern  int harmonic_compute(double *serie, double *residual, int n, harmonic_t x, spectrum_t s, int maxstep, double *a, double *G);
extern void harmonic_compute(double *serie, double mask,  double *residual, int n, harmonic_t x, spectrum_t s);


//from woce.c
/*
NOT DEFINED !!!
extern void woce_load(char *filename, double *h, double* t, double *mask, int *n, double *model, int nmod, char *outfile);
*/
extern int senetosa_load(const char *filename, double t1, double t2, double **h, double **t, double mask, int *n, double *dt);

extern int load_serie(const char *filename, double **time, double **z, char unit);
extern int load_serie1D(const char *, double **, double **,  char, mooring_t *);
extern int load_serie2D(const char *, double **, double ***, char, mooring_t *);

extern int matroos_load(const char *filename, double t1, double t2, double **h, double **t, double *mask, int *n, double *dt);
extern int hawaii_load(const char *filename, mooring_t *, double **elevation, double **time, double *mask, int *n, char *outfile);

extern int sample_load(const char *filename, double **h, double **u, double **v, double **p,double **t, double *mask, int *n, double *dt);
extern int sample_load(const char *filename, tseries_t & serie, double *dt);

extern int sample_save(const char *filename, double *h, double *u, double *v, double *p, double mask, double *t, int n);
extern int matroos_save(const char *filename, mooring_t mooring, double *h, double *t, double mask, int n);
extern int save_timeserie(const char *filename, double *h, double mask, double *t, int n);

extern int timeserie_save(const char *filename, mooring_t mooring, double *h, double mask, double *t, int n, int fmt);

extern int skipHeader(FILE *fic_data);
extern int aktarus_load_serie(double *time, double *sealevel, double *lon, double *lat,FILE *fic_data,int heure);
extern serie_t aktarus_reduce(serie_t & source,date_t start=NADate,date_t end=NADate);

extern double deltaOmega2separation(double deltaOmega);
extern double printAndGetSpectrumDetails(spectrum_t AnalysisList);

extern int harmonic_start(harmonic_t *harmonic,const spectrum_t AnalysisList);
extern int harmonic_start(harmonic_t *harmonic,const size_t *nndes, int nrhs);
extern int harmonic_start(harmonic_t *harmonic,int nframe);
extern harmonic_t harmonic_start(spectrum_t AnalysisList,const size_t *nndes, int nrhs, int nframe=0);

extern int harmonic_init(harmonic_t *x,double *time,int nodal,const astro_angles_t &astro_angles, double averaging_time);
extern int harmonic_init_new(harmonic_t & x,const double *time, int nframes,const spectrum_t & spectrum, int nodal, double averaging_time);

extern void harmonic_free(harmonic_t harmonic);
extern void harmonic_free_rhs(harmonic_t harmonic);

extern int harmonic_coefficients(double t, spectrum_t spectrum, float *cs, float *sn, int nodal,const astro_angles_t &astro_angles);
extern int harmonic_coefficients(double t, spectrum_t spectrum, double *cs, double *sn, int nodal,const astro_angles_t &astro_angles);
extern void harmonic_correction(double time, spectrum_t AnalysisList, hconstant_t **constants,float **buffer,const size_t *nndes, int nrhs, int nodal_corrections,const astro_angles_t &astro_angles);
extern void harmonic_storage(double t, harmonic_t harmonic, int nodal_corrections, float **buffer,const astro_angles_t &astro_angles);
extern void harmonic_correction(const harmonic_t harmonic, int frame, hconstant_t **constants, float **buffer);
extern void harmonic_storage(const harmonic_t& harmonic, int frame, float** buffer);

extern int addwave(spectrum_t *list, tidal_wave wave);
extern int removewave(spectrum_t *list, int i);
extern int spectrum_init(spectrum_t *list);
extern void spectrum_terminate(spectrum_t list);

extern void harmonic_init(spectrum_t WaveList, spectrum_t &, int nndes) __attribute__((safe_deprecated("allocates a global variable : bharm")));
extern void harmonic_save(spectrum_t &, double **b[3], int nndes, date_t start,date_t final);
extern void harmonic_end(spectrum_t &, date_t start,date_t final,int nndes, int count);
extern void harmonic_storage_obsolete(spectrum_t &, double t,int nndes,float *hmean,float *buffer[3],int *count, const astro_angles_t &astro_angles) __attribute__((safe_deprecated("Uses global variables bharm and Aharm")));

extern int decode_name(date_t actual, const char *varname, const char *name_template, char **filename);
extern int decode_name(date_t actual, const char *varname, const char *name_template, char **filename, int *);

extern void print_decode_atlasname_help(int mode=3);
extern char *decode_atlasname(const char *convention,const char *wave, const char*var, int mode=3);
extern void print_tide_decode_atlasname_help(int mode=3);
extern int tide_decode_atlasname(const char* atlas_directory, const char* atlas_convention, const char* wave, const char* varname, char** filename, int mode=3);

extern int tide_atlas2mgr(const char *filename, const char *v1, const char *v2, vector<mgr_t> mgr, int nmgr,double *a,double *G, double mask,atlas_grid_or_mesh_t *gm);
extern int tide_atlas2mgr(const char *filename, const char *v1, const char *v2, vector<mgr_t> mgr, int nmgr,float *a,float *G, float mask,atlas_grid_or_mesh_t *gm);

extern int tide_atlasSG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, float *a,float *G, float mask,const grid_t &grid,int verbose=1,const char *wave=NULL,int strict=0);
extern int tide_atlasSG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, double *a,double *G, double mask,const grid_t &grid,int verbose=1,const char *wave=NULL,int strict=0);
extern int tide_atlasSG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, complex<double> *z,complex<double> cmask,const grid_t &grid,int verbose=1,const char *wave=NULL,int strict=0);
extern int tide_atlasSG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, complex<float> *z,complex<float> cmask,const grid_t &grid,int verbose=1,const char *wave=NULL,int strict=0);

extern int tide_atlas2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, float *a,float *G, float mask,atlas_grid_or_mesh_t *gm,int verbose=1,const char *wave=NULL,int strict=0);
extern int tide_atlas2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, double *a,double *G, double mask,atlas_grid_or_mesh_t *gm,int verbose=1,const char *wave=NULL,int strict=0);
extern hconstant_t *tide_atlas2positions(const char *atlPathConv,const spectrum_t &WaveList,const char *ampNameConv,const char *phaNameConv,double *lon,double *lat,int n,const double mask,int verbose,int strict=0);
extern void scale_constants(hconstant_t *constants,int pn,int wn,double factor);
extern int tide_atlas2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, complex<double> *z,complex<double> cmask,atlas_grid_or_mesh_t *gm,int verbose=1,const char *wave=NULL,int strict=0);
extern int tide_atlas2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, complex<float> *z,complex<float> cmask,atlas_grid_or_mesh_t *gm,int verbose=1,const char *wave=NULL,int strict=0);

extern hconstant_t atlas_constants(const char *atlas_convention,const spectrum_t &s, double x, double y);

extern hconstant_t * load_atlas(const string & atlas_template,const string & v1,const string & v2, const spectrum_t & WaveList, int *pn_=0, poc_var_t *bv=0);
extern int load_atlas(const char *filename,const char *v1,const char *v2, grid_t *grid, complex<float> **tide, complex<float> *cmask, int mode, int verbose);
extern int swap_XY2YX(grid_t grid, complex<float> *buffer);
extern int swap_XY2YX(grid_t grid, float *buffer);
extern int swap_XY2YX(grid_t grid, double *buffer);

extern int eot_loadc1(const char* filename,grid_t *grid, complex<float> **tide, complex<float> *mask);
extern int sirocco_loadc1(const char* filename,grid_t *grid, complex<float> **tide, complex<float> *mask);
extern int tpxo_loadc1(const char* filename,grid_t *grid,complex<float> ***tide, complex<float> *mask, int *nbuffers, spectrum_t *s);
extern int tpxo9_loadc1(const char* filename,const char* gridfile,grid_t *grid,complex<float> ***tide, complex<float> *mask, int *nbuffers, spectrum_t *s);
extern int boy_loadc1(const char* filename,grid_t *grid,complex<float> **tide, complex<float> *mask);
extern int osu_loadc1(const char* filename, grid_t *grid, complex<float> **tide, complex<float> *mask);

extern void ellipse_Madmp(const complex<double> u,const complex<double> v,double *M,double *a=NULL,double *d=NULL,double *m=NULL,double *p=NULL);
extern void ellipse_Madmp3D(const complex<double> u,const complex<double> v,const complex<double> w,double *M,double *a=NULL);
extern ellipse_t ellipse_parameter(const complex<double> vx, const complex<double> vy);
extern ellipse_t ellipse_parameter(const complex<float> vx, const complex<float> vy);
extern void ellipse_parameter(complex<float> vx,complex<float> vy,float *rmin,float *rmax,float *pol, float *dir, float *time);

extern int tides_savec1(const char *output, grid_t grid, complex<float> *tide,complex<float> cmask, float scale);
extern int tides_savec1(const char *output, grid_t grid, fcomplex *tide, fcomplex cmask, float scale, string & Aname, string & Gname);
extern int tides_savec1(const char *output, const char *v1, const char *v2, const char *name, const char *units, grid_t zgrid, complex<float> *cbuf, complex<float> cmask, pocgrd_t ncgrid);
extern int tides_savec1(const char *output, const char *v1, const char *v2, const char *name, const char *units, grid_t grid, complex<float> *cbuf, complex<float> cmask);

extern tidal_wave wave_from_code(int code);
extern double omega_from_frequency_or_period(const float fOrT, const char *unit);
extern tidal_wave wave_from_doodson(const int doodson);
extern int doodson_from_wave(const tidal_wave & w);

extern double tide_omega(tidal_wave wave, const char *units="dph");
extern double tide_period(tidal_wave wave, bool initIf0=true);

extern void alias_frequency(double alias, double *omega);

//from tides-friction.cpp

#define TidesHarmonicConstituentsModeling_REF "Provost, Rougier and Poncet 1981 http://dx.doi.org/10.1175/1520-0485(1981)011%3C1123%3ANMOTHC%3E2.0.CO%3B2"

extern void frodom(complex<double> v1,complex<double> v2,double *ro,double *roprim);
extern void frosec(complex<double> v1,complex<double> v2,double *r,double *rprim,double *rsecon);

extern void spectral_friction01(zmatrix2x2_t *fric,complex<double> u,complex<double> v);
extern void spectral_friction02(zmatrix2x2_t *fric,complex<double> u,complex<double> v);

extern int harmonic_check(spectrum_t, double, double);
extern int harmonic_check(spectrum_t, double, double, double *);

extern int harmonic_optimize(spectrum_t, double, double, double *);

#include "admittance-mts.h"

#endif
