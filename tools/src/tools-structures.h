
/*******************************************************************************

  T-UGO tools, 2006-2018

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file
*/
/*----------------------------------------------------------------------------*/

#if TOOLS_STRUCTURES_H == 0
#define TOOLS_STRUCTURES_H 1

#include <netcdf.h>             //for nc_type,NC_...
#include <string.h>             //for str...
#include <cfloat>               //for DBL_MAX
#include <limits>               //for numeric_limits
#include <iostream>             //for cout
#include <algorithm>

#define BinaryRead(vbl,type,n,file) bread_dg((char *)(&(vbl)),(sizeof(type)),(n),(file))
#define BinaryWrite(vbl,type,n,file) bwrite_dg((char *)(&(vbl)),(sizeof(type)),(n),(file))
#define NINT(a)   (floor(a+0.5))
#define nint(a)   (int) (floor(a+0.5))
#define ST(a,b,i) ((a+(i)-1)->b)

#include <solverlib.h>

extern char *poc_strdup(const char *source);

#include "poc-array.hpp" //for deletep()

#undef PREFIX
#undef INITIALVALUE
#ifdef MAIN_SOURCE
#   define PREFIX
#   define INITIALVALUE(x) =x
#else
#   define PREFIX extern
#   define INITIALVALUE(x)
#endif

/*----------------------------------------------------------------------------
 parallel computing variables (MPI) */
PREFIX int gParallel INITIALVALUE(0);
PREFIX int gCPU_ID INITIALVALUE(0), gCPU_MASTER INITIALVALUE(0), gCPU_STREAMER INITIALVALUE(0);
PREFIX int gnCPUs INITIALVALUE(1);

/* NOTE: DO NOT ENABLE CODE BELOW AS IT WILL BREAK POCViP !!!
#undef PREFIX
#undef INITIALVALUE
*/


typedef struct {
  double **u,  **v, **w;         /* velocities                       */
  double **U,  **V, **W;         /* massic velocities  (u x rho) ?   */
  double **T,  **S, **rho, **p;  /* tracers, density and pressure    */
  double **Ex, **Ey;             /* Horizontal diffusion             */
  double **Ax, **Ay;             /* Horizontal advection             */
  double **Kh3D;                 /* Horizontal diffusion coefficient */
  double **px, **py   ;          /* H-baroclinic pressure gradient   */
  double **z;                    /* Immersion                        */
  double *elevation;             /* sea surface elevation            */
  double *umean,*vmean;          /* mean 2D velocities over a 3D step*/
  int    umode,wmode,tmode;
} ugo_state_t;

// struct List {
//   void *data;
//   List *next;
//   List *prev;
// };

/*--------------------------------------------------------------------------
  Node based arrays are part of structure type "node"
  model2D.variables at time step k+1 (Hp etc), k (H etc) and k-1 (Hm etc)
---------------------------------------------------------------------------*/

typedef struct {
  double h,dhdx,dhdy,r;     /* mean depth, gradient vector and rugosity   */
  double Hp,H,Hm;           /* next, current, previous  depth             */
  float  zb,zbm;            /* bottom motion                              */
  float  up,vp,u,v,um,vm;   /* next, current, previous  u and v currents  */
  double hmean,umean,vmean,humean,hvmean; /* mean h, current and transpor */
  double lon,lat;           /* longitude, latitude                        */
  float  tau,Ex,Ey;         /* friction & viscosity coefficients          */
  float  ux,uy,vx,vy;       /* velocity gradient                          */
  float  vic;               /* variable ice cover coefficient             */
  float  vfc,vbv;           /* variable friction coefficient, background velocity */
  float  lmm;               /* lumped mass matrix                         */
  int    sponge,stable;     /* a spongious barrier flag, stability flag   */
  float  C,S;               /* cosine and sine of the latitude            */
  float  N;                 /* Mean Brunt-Vassala frequency               */
  float  *Ha,*Hg,*Ua,*Ug,*Va,*Vg; /*  Computed tidal constants            */
  int    code,index;        /*   node code, index in main node set        */
} state2D_t;

/*--------------------------------------------------------------------------
  code= 0 for interior nodes;
  code= n for boundary node (nth boundary, 1 is external boundary)
  code=-1 for open boundary conditions (presumedly subset of external boundary)
---------------------------------------------------------------------------*/

typedef struct {
  double Pa, Pam,tl, tlm;     /* atmospheric pressure, tidal loading      */
#ifdef PERTURB
  double dPa;                 /*   */
#endif
  complex<float>  *tlc;               /* tidal loading, true real & imaginary part*/
  float *tpA,*tpG;            /* tidal potential                          */
  float  tpx,tpy,tpxm,tpym;   /* tidal potential gradient                 */
  float  wsx,wsy,wsxm,wsym;   /* Wind stress                              */
  float  wdx,wdy,wdxm,wdym;   /* wave drag                                */
} forcing_t;

  typedef struct {
/*--------------------------------------------------------------------------
  geometry parameters*/
  double h, dhdx, dhdy;
  double lon, lat;      /* longitude, latitude                  */
  float C, S;   /* cosine and sine of the latitude      */
  float lmm;    /* lumped mass matrix                   */
/*--------------------------------------------------------------------------
  geometry parameters*/
  int stable, sponge;   /*         */
  int bcode, index, bflag;      /*         */
/*--------------------------------------------------------------------------
  dynamic parameters*/
  float vfc, r;
  double u0, chezy;
  float *Ha, *Hg, *Ua, *Ug, *Va, *Vg;   /*  Computed tidal constants       */
} parameter_obsolete_t;

#define TIDAL_WAVE_NAME_LEN 20
class tidal_wave {
private :
public :
  float Ap;                     /*  specify tidal potentiel amplitude                  */
  int order;                    /*  specify tidal potentiel development order          */
  int nT, ns, nh, np, np1; /* specify the main wave fequency V. See greenwhich_argument() */
  int nksi, nnu0, nnu1, nnu2;   /*  specify nodal argument u. See nodal_phase()        */
  int Qu, Ra;                   /*  specify nodal argument u. See nodal_phase()        */
  float shift;                  /*  See greenwhich_argument()                          */
  int formula;                 /* formula index for nodal factor f. See nodal_factor() */
  int code;                     /*  code for JMM    */
  double omega,aliased;         /*  pulsation (deg/h, see tide_omega() ) */
  char name[TIDAL_WAVE_NAME_LEN];/*                                                     */

  void reset (){
    Ap=0;
    nT=ns=nh=np=np1=0;
    nksi=nnu0=nnu1=nnu2=0;
    Qu=Ra=0;
    shift=0;
    formula=0;
    //code=0;
    omega=NAN;
    name[0]=0;
    }
  
  /* DO NOT UNCOMMENT THIS, unless you change tides.def
  tidal_wave(){
    reset();
    }
  */
  void init(const char *units="dph") {
    extern double tide_omega(tidal_wave wave, const char *units="dph");
    if(isnan(Ap)) return;
    omega = aliased = tide_omega(*this,units);
    }
};



class mgr_data_t {
  private :
  public :
  double amp, phi;
  double error;
  double bg_contamination_error;
  tidal_wave constituent;

  mgr_data_t(){
    amp   = 0.;
    phi   = 0.;
    error = 0.;
    bg_contamination_error = 0.;
    }

  mgr_data_t &operator = (const mgr_data_t &data) {
    amp   = data.amp;
    phi   = data.phi;
    error = data.error;
    bg_contamination_error = data.bg_contamination_error;
    constituent=data.constituent;
    strcpy(constituent.name, data.constituent.name);
    return(*this);
    }

};

class mgr_loc_struct {
private :
public :
  double lat,lon,immersion,depth;
  char   *units;
  
  mgr_loc_struct():
  lon(1.e+35),
  lat(1.e+35),
  immersion(1.e+35),
  depth(1.e+35),
  units(0)
  {}; // default constructor
  
  mgr_loc_struct(const mgr_loc_struct &data) {// copy constructor
    lat = data.lat;
    lon = data.lon;
    immersion = data.immersion;
    depth = data.depth;
    if(data.units!=0) units = strdup(data.units);
    else units = 0;
    }
  
  mgr_loc_struct &operator = (const mgr_loc_struct &data) { // assignment operator
    lat = data.lat;
    lon = data.lon;
    immersion = data.immersion;
    depth = data.depth;
    if(units != 0)
      free(units);
    if(data.units!=0) units = strdup(data.units);
    else units = 0;
    return(*this);
  }


};

class mgr_t {
  private :
  public :

  mgr_data_t *data;
  mgr_loc_struct  loc;
  char name[128];
  char origine[256];
  char validation[256];
  int number;
  int mindex,track;
  int duree;
  char debut[11],fin[11];
  int nwave;

  mgr_t() {
    data=0;
    number=-1;
    mindex=1;
    track=-1;
    strcpy(debut,"");
    strcpy(fin,"");
    duree=-1;
    nwave=-1;
    strcpy(validation,"??");
    loc.lon = loc.lat = loc.immersion = loc.depth = 1.e+35;
    loc.units = 0;
    }

  mgr_t(const mgr_t &mgr) {
    data  = mgr.data;
    loc   = mgr.loc;
    strcpy(name,mgr.name);
    strcpy(origine,mgr.origine);
    strcpy(validation,mgr.validation);
    number=mgr.number;
    mindex = mgr.mindex;
    track=mgr.track;
    duree=mgr.duree;
    strcpy(debut,mgr.debut);
    strcpy(fin,mgr.fin);
    nwave=mgr.nwave;
    }

  mgr_t &operator = (const mgr_t &mgr) {
   data  = mgr.data;
   loc   = mgr.loc;
   strcpy(name,mgr.name);
   strcpy(origine,mgr.origine);
   strcpy(validation,mgr.validation);
   number=mgr.number;
   mindex = mgr.mindex;
   track=mgr.track;
   duree=mgr.duree;
   strcpy(debut,mgr.debut);
   strcpy(fin,mgr.fin);
   nwave=mgr.nwave;

   return(*this);
   }
   
  int wave_index(const char *name){
    using namespace std;
    for(int i = 0; i < nwave; i++){
      if(strcmp(name, data[i].constituent.name) == 0){
        return(i);
        }
      }
//     printf("no match, try uppercase comparison\n");
    for(int i = 0; i < nwave; i++){
      string s1=name;
      string s2=data[i].constituent.name;
      transform(s1.begin(), s1.end(), s1.begin(), ::toupper);
      transform(s2.begin(), s2.end(), s2.begin(), ::toupper);
      if(s1==s2){
        return(i);
        }
      }
    return(-1);
    }
    
  void clear(){
    number = -1;
    mindex = -1;
    track = -1;
    strcpy(debut,"");
    strcpy(fin,"");
    duree = -1;
    nwave = -1;
    strcpy(validation,"??");
    deletep(&data);
    }
};


class ellipse_t {
private :
public :
/*----------------------------------------------------------------------------------------------------
  geometry parameters*/
  double
    a,///<maximum current
    b,///<minimum current
    polarisation,///<
    inclination,///<direction of maximum current in radians
    phase;///<phase of maximum current in radians

  ellipse_t() {
    a=b=polarisation=inclination=phase=0.0;
    }
};

typedef struct {
  double H, z, hu, hv;
  double pb, P;
  double u, v;
  double ux, uy, vx, vy;
  double tau, Ex, Ey, Ah;  /*        */
  float  vbv, vic, vie;    /* variable background velocity, ice cover , ice elasticity*/
  float  N;                 /* Mean Brunt-Vassala frequency         */
} state2d_obsolete_t;

class state2d_t {
private :
public :
  double *H, *z;                       /* total depth, elevation, transport                             */
  double *H__,*z__;                    /* total depth at momentum nodes                                 */
  double *u__,*v__;                    /* velocities at elevation nodes                                 */
  double *U__,*V__;                    /* transport at elevation nodes                                  */
  double *pb, *P;                      /* bottom pressure, integrated presssure                         */
  double *u, *v,  **w, *wres;          /* velocities                                                    */
  double *U, *V;                       /* transport                                                     */
  double *ux, *uy, *vx, *vy;           /* velocities gradients                                          */
  int    *suspect, *stability;         /* stability control                                             */
  double *rtau,*rx, *ry, *taux, *tauy; /* friction                                                      */
  double *taux__, *tauy__;             /* friction at elevation nodes                                   */
  double *Ex, *Ey, *Ah, *Kh, *lmts;    /* diffusion, diffusion coefficient, lateral mixing time scale   */
  double *Ex__, *Ey__, *Ah__;          /* diffusion, diffusion coefficient, lateral mixing time scale   */
  double *Ax, *Ay;                     /* diffusion, diffusion coefficient, lateral mixing time scale   */
  double *delta3Dx, *delta3Dy;         /* mean 3D contribution corrections                              */
  float  *vic, *vie;                   /* ice cover, ice elasticity                                     */
  float  *N;                           /* Mean Brunt-Vaisala frequency                                  */
  int    colocation;                   /* flag indicating pressure and velocity colocation              */

  state2d_t() {
    H=0; z=0;
    H__=0; z__=0;
    pb=0;P=0;
    u=0; v=0; w=0; wres=0;
    U=0; V=0;
    u__=0; v__=0;
    U__=0; V__=0;
    ux=0; uy=0; vx=0; vy=0;
    suspect=0; stability=0;
    rtau=0; rx=0; ry=0; taux=0; tauy=0;
    Ex=0; Ey=0; Ah=0; Kh=0; lmts=0;
    Ax=0; Ay=0;
    delta3Dx=0; delta3Dy=0;
    vic=0;vie=0;
    N=0;
    Ex__=0; Ey__=0; Ah__=0;
    taux__=0; tauy__=0;
    colocation=-1;
    }

  void import_z(int target, int source, state2d_t state2D) {
    H[target]=state2D.H[source];
    z[target]=state2D.z[source];
    P[target]=state2D.P[source];
    pb[target]=state2D.pb[source];
    }

  void import_u(int target, int source, state2d_t state2D) {
    u[target]=state2D.u[source];
    v[target]=state2D.v[source];
    V[target]=state2D.U[source];
    v[target]=state2D.V[source];
//     ux[target]=state2D.ux[source];
//     uy[target]=state2D.uy[source];
//     vx[target]=state2D.vx[source];
//     vy[target]=state2D.vy[source];
    Ex[target]=state2D.Ex[source];
    Ey[target]=state2D.Ey[source];
/**----------------------------------------------------------------------------
  Ah is implicitely u node difusion coefficient, Kh is u-gradient node value */
    Ah[target]=state2D.Ah[source];
    Ax[target]=state2D.Ax[source];
    Ay[target]=state2D.Ay[source];
    delta3Dx[target]=state2D.delta3Dx[source];
    delta3Dy[target]=state2D.delta3Dy[source];
    vic[target]=state2D.vic[source];
    vie[target]=state2D.vie[source];
    N[target]=state2D.N[source];
    }

  void import_g(int target, int source, state2d_t state2D) {
    ux[target]=state2D.ux[source];
    uy[target]=state2D.uy[source];
    vx[target]=state2D.vx[source];
    vy[target]=state2D.vy[source];
/**----------------------------------------------------------------------------
  Ah is implicitely u node difusion coefficient, Kh is u-gradient node value */
    Kh[target]=state2D.Kh[source];
    }

  void destroy() {
    deletep(&H);
    deletep(&z);
    deletep(&pb);
    deletep(&P);
    deletep(&u);
    deletep(&v);
    deletep(&U);
    deletep(&V);
    deletep(&ux);
    deletep(&uy);
    deletep(&vx);
    deletep(&vy);
    deletep(&Ex);
    deletep(&Ey);
    deletep(&rtau);
    deletep(&rx);
    deletep(&ry);
    deletep(&taux);
    deletep(&tauy);
    deletep(&Ah);
    deletep(&Kh);
    deletep(&lmts);
    deletep(&Ax);
    deletep(&Ay);
    deletep(&delta3Dx);
    deletep(&delta3Dy);
    deletep(&vic);
    deletep(&vie);
    deletep(&N);
    switch(colocation) {
      case 0:
        deletep(&H__);
        deletep(&z__);
        deletep(&u__);
        deletep(&v__);
        deletep(&U__);
        deletep(&V__);
        deletep(&Ex__);
        deletep(&Ey__);
        deletep(&Ah__);
        deletep(&taux__);
        deletep(&tauy__);
        break;

      case 1:
        break;

      default:
        STDOUT_BASE_LINE("colocation flag %d not valid\n",colocation);
        exit(-1);
        break;
     }
   }
};

typedef struct {
  double *u, *v, *w;
  double *U, *V, *W;
  double *T, *S, *rho, *p;
  double *Ex, *Ey;      /* Horizontal diffusion         */
  double *Ax, *Ay;      /* Horizontal advection         */
  double *N;
  double *z;
  double H, tau;
} state3d_t;

typedef struct {
  double Pa;            /* atmospheric pressure                     */
  double tl, zb;        /* tidal loading,bottom motion              */
#ifdef PERTURB
  double dPa;           /* atmospheric pressure perturbation        */
#endif
  double prx, pry;      /* pressure gradient                         */
  double tpx, tpy;      /* tidal potential gradient                  */
  double wsx, wsy;      /* Wind stress                               */
  float  pax, pay;      /* atmospheric pressure gradient             */
  double wdx, wdy;      /* internal wave drag                        */
  double Sxx, Sxy, Syy; /* surface (ocean) wave drag                 */
//  fcomplex *tlc;      /* tidal loading, true real & imaginary part */
  float *tpA, *tpG;     /* tidal potential                           */
} action_obsolete_t;

class action_t {
private :
public :
  double *TFlux;                 /* Temperature flux                          */
  double *SFlux;                 /* Salinity flux                             */
  double *Pa;                    /* atmospheric pressure                      */
  double *LSA, *TidalPotential;  /* tidal loading, potential                  */
  double *BottomDeformation;     /* bottom motion                             */
  double *prx, *pry;             /* pressure gradient                         */
  double *tpx, *tpy;             /* tidal potential gradient                  */
  double *wsx, *wsy;             /* Wind stress                               */
  double *wdx, *wdy;             /* internal wave drag                        */
  double *Sxx, *Sxy, *Syy;       /* surface (ocean) wave drag                 */
  double *RSdvg_x, *RSdvg_y;
#ifdef PERTURB
  double dPa;                    /* atmospheric pressure perturbation         */
#endif

  action_t() {
    TFlux=0;
    SFlux=0;
    Pa=0;
    LSA=0;
    prx=pry=0;
    tpx=tpy=0;
    wdx=wdy=0;
    wsx=wsy=0;
    Sxx=Sxy=Syy=0;
    RSdvg_x=RSdvg_y=0;
    }

  void allocate(int zdim, int udim, int tdim) {
    TFlux=0;
    SFlux=0;
    Pa=0;
    LSA=0;
    prx=pry=0;
    tpx=tpy=0;
    wdx=wdy=0;
    wsx=wsy=0;
    RSdvg_x=RSdvg_y=0;
    if(tdim!=0) {
      TFlux=new double[tdim];
      SFlux=new double[tdim];
      }
    if(zdim!=0) {
      Pa =new double[zdim];
      LSA=new double[zdim];
      }
    if(udim!=0) {
      prx=new double[udim];
      pry=new double[udim];
      tpx=new double[udim];
      tpy=new double[udim];
      wdx=new double[udim];
      wdy=new double[udim];
      wsx=new double[udim];
      wsy=new double[udim];
      RSdvg_x=new double[udim];
      RSdvg_y=new double[udim];
      }
    }

  void destroy() {
    deletep(&TFlux);
    deletep(&SFlux);
    deletep(&Pa);
    deletep(&LSA);
    deletep(&prx);
    deletep(&pry);
    deletep(&tpx);
    deletep(&tpy);
    deletep(&wdx);
    deletep(&wdy);
    deletep(&wsx);
    deletep(&wsy);
    deletep(&Sxx);
    deletep(&Sxy);
    deletep(&Syy);
    deletep(&RSdvg_x);
    deletep(&RSdvg_y);


    }
};


/* *----------------------------------------------------------------------------
  Unstructured grid */

// #include "fe-classes.h"


class meteo_t {
private :
public :
  float *P, *u, *v;            /* atmospheric pressure, wind, temperature   */
  float *rge, *drn;            /* wind range and direction                  */
  float *wsx, *wsy;            /* wind stress                               */
  float *T, *q;                /* temperature, humidity                     */
  float *Hsw, *Hlw;            /* descending heat flux short and long waves */
  float *Fw;                   /* rainfall                                  */
  float *dP;                   /* data assimilation stuff                   */

  meteo_t () {
    P=u=v=0;
    rge=drn=0;
    wsx=wsy=0;
    T=q=0;
    Hsw=Hlw=0;
    Fw=0;
    dP=0;
    }
    
  int allocate(int zdim, int udim) {
    P   = new float[zdim];
    if(P==0) return(-1);
    rge = new float[udim];
    if(rge==0) return(-1);
    drn = new float[udim];
    if(drn==0) return(-1);
    return(0);
    }
    
  void destroy() {
    deletep(&P);
    deletep(&rge);
    deletep(&drn);
    }
};


// typedef struct
// {
//   state2D_t    nodes;            /*    */
//   forcing_t forcing;          /*    */
//   int       nn;               /*    */
//   triangle_t elements;         /*    */
//   int       ne;               /*    */
// } model_t;

typedef struct {                 /*  Open boundary nodes                     */
  int node;                      /*  Node index in the neighbour list        */
  float *h_a,*h_G;               /*  Specified h tidal constants	     */
  float *u_a,*u_G;               /*  Specified u tidal constants	     */
  float *v_a,*v_G;               /*  Specified v tidal constants	     */
  float h,z,u,v;                 /*  Model solution                          */
  float h_ext,z_ext,u_ext,v_ext; /*  Specified external solution             */
  float nx,ny,size;              /*  normal vector, segment length           */
} obc2D_t;

typedef struct
{
  double tl;                  /* tidal loading                      */
  float *tlA,*tlG;            /* tidal loading                      */
  complex<float> *tlc;        /* tidal loading                      */
} loading_t;


typedef struct {              /*  Open boundary nodes                                  */
  int node;                   /*  Node index in the neighbour list                     */
  int indx;                   /*  boundary indes in upper level mesh boundary set      */
  float *a,*G;                /*  Specified tidal constants                            */
  double z,u,v,Hu,Hv,H,ib;    /*  Specified elevation and currents from external model */
  double zm,um,vm,Hum,Hvm,Hm,ibm;    /*  Specified elevation and currents from external model */
  float sine,cosine,size;     /*  Ny, Nx of normal vector, segment length              */
} z_spec;

typedef struct {              /*  land boundary nodes                     */
  int node;                   /*  node index in the neighbour list        */
  float sine,cosine,size;     /*  Ny, Nx of normal vector, segment length */
} l_bndry;

typedef struct {
  char var;
  float lon,lat;
  int ind;
} plot_points;

typedef struct {
  float lon,lat;             /*  location of the data   */
  float *a,*G;               /*  specified tidal data   */
} tide_def;


class spectrum_t {
 public:
  tidal_wave *waves;     /*  array of tidal waves              */
  int *prescribed;      /*  data: 0 from model, 1 from admittance, 2 from equilibrium */
  int n;                /*  number of waves in array          */
  int nmax;             /*  maximum number of waves in array  */
  
 private:
  /* all constructors MUST call this */
  void NULLPointers(){
    waves=NULL;
    prescribed=NULL;
    }
  
  /* this is a constructor */
  friend void init_spectrum(spectrum_t *spectrum, int *n);
  
 public:
  spectrum_t(){
    NULLPointers();
    n=0;
    nmax=0;
    }

  void init(int nn){
    n=0;
    nmax=nn+1;//to enable adding Z0 when necessary
    deletep(&waves);
    waves=new tidal_wave[nmax];
    deletep(&prescribed);
    prescribed = new int[nmax];
    }

  spectrum_t(int nn){
    NULLPointers();
    init(nn);
    }

  spectrum_t(spectrum_t reference, mgr_t mgr){
    const int nn=mgr.nwave;
    NULLPointers();
    init(nn);
    for(int k=0;k<nn;k++) {
      this->add(reference, mgr.data[k].constituent.name, 0);
      }
    }

  spectrum_t &operator = (const spectrum_t &source) {
    n = source.n;
    if(source.nmax<=0) {
      nmax = source.n;
      }
    else {
      nmax = source.nmax;
      }
    waves=new tidal_wave[nmax];
    for(size_t k = 0; k < n; k++) waves[k] = source.waves[k];
    return(*this);
    }

  void destroy() {
    deletep(&waves);
    deletep(&prescribed);
    n=0;
    nmax=0;
    }

  int wave_index(const char *name) const{
    
    if(strcasecmp(name,"NIV")==0)
      name="Z0";
    
    for(int i = 0; i < this->n; i++){
      const char *wName=waves[i].name;
      if(strcasecmp(name, wName) == 0){
        return(i);
        }      
      }
    
//     printf("no match, try uppercase comparison\n");
    for(int i = 0; i < this->n; i++){
      string s1=name;
      string s2=waves[i].name;
      transform(s1.begin(), s1.end(), s1.begin(), ::toupper);
      transform(s2.begin(), s2.end(), s2.begin(), ::toupper);
      if(s1==s2){
        return(i);
        }
      }
      
    return(-1);
    }

  int wave_index(tidal_wave w){
    if(w.name[0]==0) return(-1);
    
    return wave_index(w.name);
    }

  tidal_wave wave(const char *name){
    int i = wave_index(name);
    
    if(i>=0)
      return(waves[i]);
    
    printf("no match for %s\n", name);
    tidal_wave nullwave;
    nullwave.reset();
    
    return(nullwave);
    }

  int add(tidal_wave w,int doInit=1){
    int k;
    k=this->wave_index(w.name);
    if(k!=-1) return(k);
    if(n<nmax) {
      waves[n]=w;
//      prescribed[n]=0;
      if(doInit)
        waves[n].init();
      n++;
      return(n-1);
      }
    return(-1);
    }

  int remove(tidal_wave w){
    int k,l;
    k=this->wave_index(w.name);
    if(k==-1) return(0);
    for(l=k+1;l<n;l++) {
      waves[l-1]=waves[l];
      if(prescribed!=0) prescribed[l-1]=prescribed[l];
      }
    n--;
    return(0);
    }

  int remove(int k){
    int l;
    if(k<0) return(0);
    if(k>n-1) return(0);
    for(l=k+1;l<n;l++) {
      waves[l-1]=waves[l];
      if(prescribed!=0) prescribed[l-1]=prescribed[l];
      }
    n--;
    return(0);
    }

  int add(spectrum_t reference, const char *name, int allowFile){
  ///add a wave found in a list
  /**
  \param reference list of waves
  \param *name wave name
  \returns the index of the wave or -1 if error
  */
extern tidal_wave wave_from_doodson(const int doodson);
extern double omega_from_frequency_or_period(const float fOrT, const char *unit);
extern spectrum_t spectrum_init_ref(string refClass, bool Z0, bool exitOnError=true);

    int k;
    tidal_wave w=reference.wave(name);
    
    if(w.name[0]==0){///If the wave is not found in the list:
      int l=strlen(name);
      char *endName;
      do{
        ///-it checks if it is a Doodson number
        if(l<6 || 8<l)break;
        int doodson=strtol(name,&endName,10);
        if(*endName!=0 || doodson<0)break;//invalid characters in name
        ///and if so, it gets the wave by calling wave_from_doodson()
        w=wave_from_doodson(doodson);
        }while(0);
      
      do{
        if(w.name[0]!=0)break;//it was a Doodson number
        
        ///-or it checks if it is a number
        float fOrT=strtof(name,&endName);
        
        if(endName==name){//invalid characters in name
          
          if(not allowFile) return -1;

          bool Z0=false;
          spectrum_t ref=spectrum_init_ref(name,Z0,false);
          
          if(ref.n==0)
            return -1;
          
          duplicate(ref);
          return n;
          }
        
        ///and if so, it gets the wave by calling wave_from_frequency_or_period().
        w.omega=omega_from_frequency_or_period(fOrT,endName);
        strcpy(w.name,name);
        w.Ap=NAN;//de-activate initialisation
        }while(0);
      }
    
    k=this->wave_index(name);
    ///If the wave is already in the list, it returns its index.
    if(k!=-1) return(k);
    ///If there is space for it, it adds the wave found in the list.
    if(n<nmax) {
      k=n;
      n++;
      waves[k]=w;
//      prescribed[n]=0;
      waves[k].init();
      return k;
      }
    return(-1);
    }

  void init(spectrum_t reference, char **wavename, int nw=0){
  ///initialises a list of waves
  /**
  \param reference list of all known waves
  \param **wavename NULL-terminated array of wave names
  \param nw if wavename is not NULL-terminated, specify the number of waves as 3rd argument.
    If \c nw is negative, add \c -nw to the allocated number of waves.
  */
    int j,extra=0;
    if(nw<=0){
      extra=-nw;
      for(nw=0;wavename[nw]!=NULL;nw++);
      }
    init(nw+extra);
    for(j=0;j<nw;j++) {
      if(add(reference,wavename[j],1)==-1)
        fprintf(stderr,"wave %s does not exist\n",wavename[j]);
      }
    }

  void init(const spectrum_t &reference, char *wavename){
  ///initialises a wave
  /**
  \param reference list of all known waves
  \param *wavename wave name
  wrapper for init(spectrum_t,char**,int)
  */
    init(reference,&wavename,1);
    }

  void duplicate(spectrum_t reference){
    int j;
    init(reference.n);
    for(j=0;j<reference.n;j++) {
      add(reference.waves[j]);
      }
    }

  void waves_init(){
  /// Ensures all waves[].omega are initialised
    for(int i = 0; i < n; i++) {
      if(waves[i].omega == 0.) {
        waves[i].init();
        }
      }
    }
};




class data_t
{
 private:
  size_t nvalues;

 public:
  double   time;
  double   lon,lat;
  int      cycle, pos;
  float    values[30];
  int      used;

  data_t():used(0){for(size_t k = 0; k != 30; ++k){values[k] = static_cast<float>(0);}};// class constructor
  ~data_t(){used = 0;}; // class destructor
  data_t(const data_t& source){ // copy constructor
    lon = source.lon;
    lat = source.lat;
    time = source.time;
    cycle = source.cycle;
    used = source.used;
    for(size_t k = 0; k < 30; ++k) values[k] = source.values[k];
  };

  data_t &operator = (const data_t &data) { // assignment operator
    lon = data.lon;
    lat = data.lat;
    time = data.time;
    cycle = data.cycle;
    used = data.used;
    for(size_t k = 0; k < 30; k++) values[k] = data.values[k];
    return(*this);
  }


};

class serie_t {
public:
  data_t *data;
  int    count;
  float  mask;
  
  serie_t():count(0), mask(0), data(0) {};
  
  ~serie_t(){if(data != 0) {delete [] data;}};
  
  serie_t(const serie_t& source){
    count = source.count;
    mask = source.mask;
    data = new data_t[source.count];
    for(size_t i = 0; i < source.count; ++i){
      data[i] = source.data[i];
    };
  }

  serie_t& operator= (const serie_t& source){ // assignment operator
  count = source.count;
  mask = source.mask;
  data = new data_t[source.count];
  for(size_t i = 0; i < source.count; i++){
    data[i] = source.data[i];
    }
  return(*this);	
  }
 
};
 

#define null_date date_t (0,0,0,0.)


typedef struct {
  char    comment[4][80];
  int     ni,nj,nk,nt,nd,nv,code;
  float   xmin, ymin,dx, dy;
  float   spec;
  float  *levels;
  size_t  size;
  } bmgheader_t;

typedef struct {
  char    name[60];
  double  t,p;///<longitude and latitude
  float   h;
  int     code,nwave;
  char    wave[100][10];
  char    filename[256],time_zone[32];
  complex<float> elevation[100];
  } pressure_station;

class mooring_t {
public:
  char    *name;
  int     type;
  int     code;
  double  lon,lat;
  float   depth, immersion;
  bool    initialized;
  
  mooring_t() {
    name=0;
    type=-1;
    code=-1;
    lon=lat=1.e+35;
    depth=-1;
    initialized=false;
    }
    
  mooring_t(mgr_t mgr) {
    this->name=strdup(mgr.name);
    lon=mgr.loc.lon;
    lat=mgr.loc.lat;
//    immersion=mgr.loc.immersion;
    depth=mgr.loc.depth;
    initialized=true;
    }
  };

typedef struct {
  char z,u,v;
  } archive1_t;

typedef struct {
  short z,u,v;
  } archive2_t;

typedef struct {
  float z,u,v;
  } archive4_t;

typedef struct {
  double z,u,v;
  } archive8_t;


typedef struct {
  double t,p;
  float  h;
  int    elt,node[6];
  double beta[6];
  } brozouf;

  
/*--------------------------------------------------------------------------
  basic_t was previously define in the define of archive
---------------------------------------------------------------------------*/

typedef struct {
  int    element;
  int    nodes[3];
  double beta[3];
  } basic_obsolete_t;


#undef PREFIX

#endif
