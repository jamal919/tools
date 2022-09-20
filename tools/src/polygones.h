
#ifndef POLYGONES_H
#define POLYGONES_H

#include "functions.h"
#include "map-classes.h"
#include "poc-netcdf.hpp"
#include "datastream.h"
#include "maths.h"

#define PLG_LINES_PARALLEL               -2
#define PLG_LINES_NOT_SECANT             -1
#define PLG_LINES_SECANT                  0
#define PLG_LINES_SECANT_AT_EXTRIMITY     1
#define PLG_LINES_IDENTICAL               2

#define PLG_FORMAT_UNKNOWN     -1
#define PLG_FORMAT_BIN          0
#define PLG_FORMAT_ASCII        1
#define PLG_FORMAT_SCAN         2
#define PLG_FORMAT_GSHHS        3
#define PLG_FORMAT_MIF          4
#define PLG_FORMAT_VCT00        5
#define PLG_FORMAT_XISO         6
#define PLG_FORMAT_BOUNDARIES   7
#define PLG_FORMAT_SHP          8
#define PLG_FORMAT_NETCDF       9
#define PLG_FORMAT_XY          10
#define PLG_FORMAT_MAPINFO     11
#define PLG_FORMAT_HISTOLITT   12
#define PLG_FORMAT_KML         13
#define PLG_FORMAT_POCFORMULA  14

#define PLG_FORMAT_GSHHS_1_6    0
#define PLG_FORMAT_GSHHS_2_0    1

//DEFINITION DE CONSTANTES
//int plg_point_interior=1,plg_point_boundary=0,plg_point_exterior=-1;
#define PLG_POINT_EXTERIOR -1
#define PLG_POINT_BOUNDARY  0
#define PLG_POINT_INTERIOR  1


/*------------------------------------------------------------------------------
  KEEP BOTH OF THESE TRUE :
    PLG_INIT_TP == PLG_SPHERICAL
    PLG_INIT_XY == PLG_CARTESIAN */

#define PLG_SPHERICAL    0
#define PLG_CARTESIAN    1
#define PLG_LOXODROMIC   2
#define PLG_ORTHODROMIC  3

enum plg_init_mode {PLG_INIT_TP,PLG_INIT_XY,PLG_INIT_SEPARATE,PLG_INIT_SHARED};
/*----------------------------------------------------------------------------*/

class plg_t;
extern int plg_single(const plg_t & polygones, double lon, double lat, int *in, int mode=PLG_CARTESIAN, int check=0);

extern int plg_isclosed(const plg_t & polygon);
extern int close_plg(plg_t *polygon);

class plg_filter_t {
private:
public:
  double dx, lengthscale, factor;
  double isocontour;
  SGfield_t<float> resolution;
  bool subdivide;
  bool enforcement, tuned;
  bool save_map, save_mesh;
  
  plg_filter_t() {
    dx=NAN;
    lengthscale=NAN;
    factor=NAN;
    isocontour=NAN;
    enforcement=false;
    tuned=false;
    save_map=false;
    save_mesh=false;
    }
    
  void init(double _scale, const SGfield_t<float> & _resolution, double _factor, double _dx, bool _subdivide, bool _tuned, bool _enforcement, bool _save_map, bool _save_mesh) {
    dx          =_dx;
    lengthscale =_scale;
    factor      =_factor;
    resolution  =_resolution;
    isocontour  =NAN;
    subdivide   =_subdivide;
    enforcement =_enforcement;
    tuned       =_tuned;
    save_map    =_save_map;
    save_mesh   =_save_mesh;
    }
  
  void check() {
/*------------------------------------------------------------------------------
    parsing : dx is cartesian working grid resolution*/
    if(dx==0.0 && resolution.x!=0) {
      double rmin;
      int m,i,j;
      rmin=amin(resolution.x,resolution.grid->Hsize(),&m);
      dx=rmin/4.0;
      resolution.grid->ij(m,&i,&j);
      printf("%s : default working grid dx not given, automatically set dx=[%d;%d]%lfm \n",__func__,i,j,dx);
      }
    
/*------------------------------------------------------------------------------
    parsing : lengthscale is smoothing typical length*/
    if(lengthscale==0) {
      lengthscale=dx*2;
      printf("%s : default smoothing lengthscale not given, automatically set lengthscale=%lf m\n",__func__,lengthscale);
      }
    }
    
  void destroy() {
    }
};


class plg_t {
private :
public :
  double *x,*y,*z;        // x and y cartesian coordinates
  double *t,*p;           // theta and phi spherical coordinates, longitude and latitude
  int    id;
  int    npt;
//  int    closed;
//   void   *data;
  char   *flag;
  
 private:
  /* all constructors MUST call this */
  void NULLPointers(){
    x=0;
    y=0;
    z=0;
    t=0;
    p=0;
//     data=0;
    flag=NULL;
    npt=0;
    id=-1;
//    closed=-1;
    }
  
 public:
  
  void setModeFromShared(int mode){
    
    if(x!=t or y!=p or t==0 or p==0)
      TRAP_ERR_EXIT(ENOEXEC,"%s called when t=%p; x=%p; p=%p; y=%p;\n", __func__, t,x,p,y);
    
    switch (mode) {
      case PLG_INIT_TP:
        x=NULL;
        y=NULL;
        break;
      case PLG_INIT_XY:
        t=NULL;
        p=NULL;
        break;
      case PLG_INIT_SEPARATE:
        x=new double[npt];
        y=new double[npt];
        valcpy(x,t,npt);
        valcpy(y,p,npt);
        break;
      case PLG_INIT_SHARED:
        break;
      default:
        TRAP_ERR_EXIT(ENOEXEC,"%s called with mode %d when only %d, %d, %d and %d are allowed\n", __func__, mode, PLG_INIT_TP, PLG_INIT_XY, PLG_INIT_SEPARATE, PLG_INIT_SHARED);
      }
    }
  
  plg_t() {
    NULLPointers();
    }
  
  plg_t(frame_t frame, int mode=PLG_SPHERICAL) {
    NULLPointers();
    init(5,PLG_INIT_SHARED);
    
    t[0]=frame.xmin;
    t[1]=frame.xmax;
    t[2]=frame.xmax;
    t[3]=frame.xmin;
    t[4]=frame.xmin;

    p[0]=frame.ymin;
    p[1]=frame.ymin;
    p[2]=frame.ymax;
    p[3]=frame.ymax;
    p[4]=frame.ymin;
    
    setModeFromShared(mode);
    }
  
  plg_t(double x1,double y1,double x2,double y2, int n=2, int mode=PLG_INIT_SEPARATE) {
    int k;
    
    NULLPointers();
    init(n,PLG_INIT_SHARED);
    
    double dx=(x2-x1)/(double) (n-1),dy=(y2-y1)/(double) (n-1);
    
    t[0]=x1;
    t[n-1]=x2;
    p[0]=y1;
    p[n-1]=y2;
    
    for(k=1;k<n-1;k++) {
      t[k]=t[0]+k*dx;
      p[k]=p[0]+k*dy;
      }
    
    setModeFromShared(mode);
    }
  
  void destroy() {
    if(x!=0)
      delete[]x;
    if(x==t){
      t=NULL;
      }
    else{
      deletep(&t);
      t=NULL;
      }
    x=NULL;
    
    if(y!=0)
      delete[]y;
    if(y==p){
      p=NULL;
      }
    else{
      deletep(&p);
      p=NULL;
      }
    y=NULL;
    
    if(z!=0) {
      delete[]z;
      z=0;
      }

    deletep(&flag);
    npt=0;
    id=-1;
    }
  
  void init(int n, const plg_init_mode shared=PLG_INIT_SHARED);
  
  void init(const plg_init_mode shared) {
    init(npt,shared);
    }
  
  bool closed() {
    int status=plg_isclosed(*this);
    return(status==1);
    }
  
  int SetFlag(char value) {
    deletep(&flag);
    if(npt<2) return(-1);
    flag=new char[npt-1];
    for(size_t k=0;k<npt-1;k++) flag[k]=value;
    return(0);
    }
  
  void duplicate(const plg_t & src) {
    this->destroy();
    this->init(0, (src.x==src.t && src.y==src.p)?PLG_INIT_SHARED:PLG_INIT_SEPARATE);
    npt=src.npt;
    poc_copy(x,src.x,npt);
    poc_copy(y,src.y,npt);
    if(src.x!=src.t && src.y!=src.p) {
      poc_copy(t,src.t,npt);
      poc_copy(p,src.p,npt);
      }
    else {
      t=x;
      p=y;
      }
    poc_copy(flag,src.flag,npt-1);
    }
  
  void duplicate(const plg_t & src, size_t start, size_t finish) {
    this->destroy();
    if(start<finish) {
      npt=finish-start+1;
      }
    else {
      npt=src.npt+finish-start;
      }
    this->init(npt,(src.x==src.t && src.y==src.p)?PLG_INIT_SHARED:PLG_INIT_SEPARATE);
    if(src.flag!=0) this->flag=new char[npt-1];
    size_t n=0;
    if(src.t!=0 && src.p!=0)
    for(size_t k=start;k<start+npt;k++) {
      size_t l=k % src.npt;
      t[n]=src.t[l];
      p[n]=src.p[l];
      if((src.flag!=0) && (n!=npt-1)) flag[n]=src.flag[l];
      n++;
      }
    n=0;
    if(src.x!=0 && src.y!=0)
    for(size_t k=start;k<start+npt;k++) {
      size_t l=k % src.npt;
      x[n]=src.x[l];
      y[n]=src.y[l];
      n++;
      }
    }
  
  void duplicate2(const plg_t & src) {
    this->destroy();

    npt=2;

    this->init(npt,(src.x==src.t && src.y==src.p)?PLG_INIT_SHARED:PLG_INIT_SEPARATE);
    if(src.flag!=0) this->flag=new char[npt-1];
    size_t n=0;
    for(size_t k=0;k<src.npt;k+=src.npt-1) {
      size_t l=k % src.npt;
      x[n]=src.x[l];
      y[n]=src.y[l];
      t[n]=src.t[l];
      p[n]=src.p[l];
//      if(src.flag!=0) flag[n]=src.flag[l];
      if((src.flag!=0) && (n!=npt-1)) flag[n]=src.flag[l];
      n++;
      }
    }
  
  double cartesian_size(int l) const{
    if(l>npt-2) return(0.0);
    double dx=this->x[l+1]-this->x[l];
    double dy=this->y[l+1]-this->y[l];
    double L=sqrt(dx*dx+dy*dy);
    return(L);
    }
  
  double cartesian_size() const{
    double L=0;
    for(int k=0;k<npt-1;k++) L+=cartesian_size(k);
    return(L);
    }
  
  double Tarea(int k) const{
    double s=0.;
    point2D_t a=point2D_t(this->t[k],this->p[k]);
    k=k+1 % npt;
    point2D_t b=point2D_t(this->t[k],this->p[k]);
    k=k+1 % npt;
    point2D_t c=point2D_t(this->t[k],this->p[k]);
    vector2D_t u(a,b),v(a,c);
    s+=u*v;
    return(s/2.0);
    }
  
  point2D_t Tbarycenter(int k, int mode=PLG_SPHERICAL) const{
    double *xx,*yy;
    switch(mode) {
      case PLG_SPHERICAL:
        xx=this->t;
        yy=this->p;
        break;
      case PLG_CARTESIAN:
        xx=this->x;
        yy=this->y;
        break;
      }
    point2D_t a=point2D_t(xx[k],yy[k]);
    k=k+1 % npt;
    point2D_t b=point2D_t(xx[k],yy[k]);
    k=k+1 % npt;
    point2D_t c=point2D_t(xx[k],yy[k]);
    vector2D_t u(a,b),v(a,c);
    point2D_t p=a+u/3.+v/3.;
    p=point2D_t((a.x+b.x+c.x)/3.0,(a.y+b.y+c.y)/3.0);
    return(p);
    }
  
  point2D_t barycenter() const{
    vector2D_t u,v;
    size_t k;
    int n=npt-1;
    if(plg_isclosed(*this))
      n--;
    point2D_t p=point2D_t(this->t[0],this->p[0]);
    if(n<=1) return(p);
    for(k=1;k<n;k++) {
      p+=point2D_t(this->t[k],this->p[k]);
      }
    p/=(double) n;
    return(p);
    }
  
    
  double area() const{
    double s=0.;
    vector2D_t u,v;
    size_t k,kk;
    double ss;
    int n=npt-1;
    if(plg_isclosed(*this))
      n--;
    point2D_t a=point2D_t(this->t[0],this->p[0]);
    for(k=1;k<n;k++) {
      point2D_t b=point2D_t(this->t[k],this->p[k]);
      kk=k+1;
      point2D_t c=point2D_t(this->t[kk],this->p[kk]);
      u=vector2D_t(a,b);
      v=vector2D_t(a,c);
      ss=u*v;
      s+=ss;
      }
    return(s/2.0);
    }
  
  double angle(int k) const{
/*------------------------------------------------------------------------------
    return consecutive segment angle (rad); aligned segments if angle = 0.0 */
    double angle=0.;
    int km, kp;
//     bool f=this->closed();
    bool f=true;
    if( (k==0) && f) {
      km=this->npt-2;
      kp=1;
      }
    else if( (k==this->npt-1) && f) {
      km=this->npt-2;
      kp=1;
      }
    else {
      km=(k-1+this->npt) % this->npt;
      kp=(k+1) % this->npt;
      }
    const vector2D_t a(this->x[k]-this->x[km],this->y[k]-this->y[km]);
    const vector2D_t b(this->x[kp]-this->x[k],this->y[kp]-this->y[k]);
    angle=a^b;
    return(angle);
    }
  
  point2D_t InsidePoint(int mode=PLG_SPHERICAL) const{
    /// HERE!!! to be reworked
    point2D_t p;
    int pos, interior;
    bool ok=false;
    double s=this->area();
    for(size_t k=0;k<npt;k++) {
      double T=this->Tarea(k);
      if((sign<double>(s)==sign<double>(T)) && fabs(T)<=fabs(s)) {
        p=Tbarycenter(k, mode);
        interior=0;
        pos=plg_single(*this, p.x, p.y, &interior, mode, 0);
        if(pos==PLG_POINT_INTERIOR) {
          ok=true;
          break;
          }
        }
      }
    if(!ok) {
      for(size_t k=0;k<npt;k++) {
        double T=this->Tarea(k);
        if((sign<double>(s)==sign<double>(T)) && fabs(T)<=fabs(s)) {
          p=Tbarycenter(k, mode);
          ok=true;
          break;
          }
        }
      }
    return(p);
    }
  
  };

class plg_array_t {
private :
public :
  plg_t *p;
  int n;

  plg_array_t() {
    p=0;
    n=0;
    }
  
  plg_array_t(size_t size) {
    n=size;
    p=new plg_t[n];
    }
  
  void destroy() {
    for(size_t k=0;k<n;k++) p[n].destroy();
    delete[] p;
    n=0;
    }
  
  };

class plg_point_t {
private :
public :
  double x,y;
  double t,p;
  int    index;

  plg_point_t() {
    x=0;
    y=0;
    t=0;
    p=0;
    index=-1;
    }
  
  plg_point_t(const plg_t & q, int k) {
    x=q.x[k];
    y=q.y[k];
    t=q.t[k];
    p=q.p[k];
    index=k;
    }
  
  void init(double x0, double y0, double t0=NAN, double p0=NAN,int index0=-1) {
    x=x0;
    y=y0;
    t=t0;
    p=p0;
    if(isnan(t))
      t=x0;
    if(isnan(p))
      p=y0;
    index=index0;
    }
  
  plg_point_t(double x0, double y0, double t0=NAN, double p0=NAN,int index0=-1) {
    init(x0,y0,t0,p0,index0);
    }
  };

typedef struct {
  int point;
  int line;
  } plg_target_t;

class line_t {
private :
public :
  plg_point_t point[2];
  line_t(double x1, double y1, double x2, double y2) {
    point[0].t=x1;
    point[0].p=y1;
    point[1].t=x2;
    point[1].p=y2;
    }
  line_t(plg_t p, int k) {
    point[0].t=p.t[k];
    point[0].p=p.p[k];
    point[1].t=p.t[k+1];
    point[1].p=p.p[k+1];
    }
  line_t(plg_t p, int k, int mode) {
    switch(mode) {
      case PLG_CARTESIAN:
        point[0].t=p.x[k];
        point[0].p=p.y[k];
        point[1].t=p.x[k+1];
        point[1].p=p.y[k+1];
        break;
      case PLG_SPHERICAL:
      default:
        point[0].t=p.t[k];
        point[0].p=p.p[k];
        point[1].t=p.t[k+1];
        point[1].p=p.p[k+1];
        break;
      }
    }
  };

class plg_interval_t {
private :
public :
  int line;
  int point[2];
  plg_interval_t() {
    line=-1;
    point[0]=-1;
    point[1]=-1;
    }
  };

class plg_block_t {
private :
public :
  plg_array_t segments;
  char  *code,type;
  int orientation;

  plg_block_t() {
    segments.n=0;
    code=0;
    }
  };

class plg_desc_t {
private :
public :
  char *scanfile;
  plg_block_t *b;
  int nblock;

  plg_desc_t() {
    b=0;
    nblock=0;
    }
  };

/*------------------------------------------------------------------------------
criteria type*/
class criteria_t {
private :
public :
  double  hmin,hmin2,hmax;             /*< topography limits (m) */
  double minratio;                     /*< topography slope limit */
  /* cell sizes are in km */
  double  minsize,maxsize,minangle;
  double  shelf_minsize,shelf_maxsize;
  double  shore_minsize,shore_maxsize;
  double  cellsize;
  double  factor1,factor2,factor3,factor4;
  double  open_shelf_limit,shelf_shore_limit;                     //depth limit between open waters and shelf
  double  limits_consistency_extent, limits_consistency_factor;
  double  maxrate;                     /*< cell density gradient limit */
  double  surface_wave_period, surface_wave_dt;
  int     mode,niterations;
  bool shelf,smooth;
  bool resample_obsolete;
  bool resample_openlimits;
  bool resample_rigidlimits;
  bool selective_open_sampling;
  vector<int> loose_openlimits;
//   bool surface_wave_cfl;
  char    descriptor   [256];
  char    boundary     [256];
  char    regulardepth [256];
  char    randomdepth  [256];
  char    delaunay     [256];
  char    rootname     [256];
  char    interior_lines [256];
  
  char    **keywords;
  
  void init() {
    hmin=hmin2=hmax=minratio=0;
    minsize=maxsize=minangle=0;
    shelf_minsize=shelf_maxsize=0;
    shore_minsize=shore_maxsize=0;
    cellsize=0;
    factor1=factor2=factor3=factor4=0;
    maxrate=0.75;
    open_shelf_limit=450;                   // default to 450m
    shelf_shore_limit=50;                   // default to 50m
    limits_consistency_extent=2.5;
    limits_consistency_factor=1.0;
    shelf_maxsize=35;
    surface_wave_period=10.0;
    surface_wave_dt=10;
    niterations=1;
    mode=0;
    shelf=smooth=0;
    resample_obsolete=resample_openlimits=resample_rigidlimits=0;
//     surface_wave_cfl=0;
    strcpy(rootname,"./");
    strcpy(regulardepth,"./depth.grd");
    descriptor[0]=boundary[0]=randomdepth[0]=delaunay[0]=0;
    interior_lines[0]=0;
    keywords=0;
    }
  
  criteria_t() {
    init();
    }
  
  int use(int n) const {
    char options=(mode >> n) % 2;
    return(options && (char) 1);
    }
    
  int uniform(double resolution) {
    init();
    cellsize=resolution;
    minsize=maxsize=minangle=resolution;
    shelf_minsize=shelf_maxsize=resolution;
    return(0);
    }
};

extern int plg_orientation (const plg_t & plg_polygones, bool debug);
extern int plg_setdirect(plg_t *polygones, int npolygones);
extern int plg_setdirect(vector<plg_t> & polygons);

extern void plg_deletep(plg_t **polygones,const int npolygones);

//from polygones-intersect.cpp
extern int plg_TestInterior(double ,double, plg_t *, int);
extern int plg_TestInterior(double, double, vector<plg_t> polygons, int mode);

extern int plg_TestInterior(const plg_t & plg_polygones,const double *lon,const double *lat, int count, char *in, bool debug);

extern char *plg_TestInterior(const double *lon,const double *lat, int count,const plg_t *polygones, int npolygones, int verbose, bool debug);
extern char *plg_TestInterior(const double *lon,const double *lat, int count,const vector<plg_t> & polygons, int verbose, bool debug);

extern char *plg_test_grid(const grid_t & grid,const plg_t *polygons, int npolygons);
extern char *plg_test_grid(const grid_t & grid, const vector<plg_t> & polygons);

extern int plg_secantpoint(double x1[2],double y1[2],double x2[2],double y2[2],double *x,double *y,double *a=0,double *b=0);
extern int plg_secantpoint (const line_t & A, const line_t & B, double *x,double *y, double *a=0, double *b=0);
extern int plg_secantpoint (const line_t & A, const line_t & B, plg_point_t & point, double *a=0, double *b=0);
extern int plg_secantpoints(const plg_t & a, const plg_t & b, point2D_t **points=0, double **angles=0, double **aindexes=0, double **bindexes=0);
extern int plg_secantpoints(const vector<plg_t> & p, const line_t line, vector<point2D_t> & points);

extern int plg_checkSecant(const vector<plg_t> & p, const vector<plg_t> & q, bool debug);

extern int plg_checkAutoSecant(plg_t *polygones, int npolygones, int mode);
extern int plg_checkAutoSecant(vector<plg_t> polygons, int mode);

extern double plg_length(plg_t plg_a, int vertex_a, int vertex_b, int mode);
extern double plg_length(plg_t plg_a, int mode);
extern double plg_minlength(plg_t plg_a, int mode);


extern int plg_find_point(plg_t *polygones, int target, int mode, double t, double p, double *d);
extern int plg_find_point(plg_t & polygon, int mode, plg_point_t & point, double *d);
extern int plg_find_point(plg_t & polygon, int mode, double t, double p, double *d);

extern int plg_NearestSegment(plg_t & polygon, double x, double y, int & m1, int & m2, double & dmin, vector2D_t & u);

extern paire_t plg_find_polygon(vector<plg_t> & polygons, plg_point_t point, double & d);
extern paire_t plg_find_polygon(vector<plg_t> & polygons, vector<int> & list, plg_point_t point, double & d);

extern int plg_NearestPolygon(plg_t *polygones, int npolygones, double x, double y, int & m1, int & m2, vector2D_t & u);
extern int plg_NearestPolygon(vector<plg_t> & polygons,         double x, double y, int & m1, int & m2, vector2D_t & u);
extern int plg_NearestPolygon(vector<plg_t> & polygons, vector<int> & exclusion, double x, double y, int & m1, int & m2, vector2D_t & u);

extern float *set_landmask01(const grid_t & grid, plg_t *polygones, int npolygones);
extern float *set_landmask01(const grid_t & grid, vector<plg_t> polygons);

extern float *set_landmask02(const grid_t & grid, plg_t *polygones, int npolygones, bool checks=false);
extern float *set_landmask02(const grid_t & grid, vector<plg_t> & polygons, bool checks);

extern float *set_landmask03(const grid_t & grid, vector<plg_t> polygons);
extern int enforce_landmask(const grid_t & grid, plg_t *polygones, int npolygones, float* mask, float value, bool tuned);

//from polygones-02.cpp
extern int plg_spherical(projPJ projection,plg_t *polygones, int npolygones);
extern int plg_spherical(projPJ ref, vector<plg_t> & polygons);

extern int plg_cartesian(projPJ projection,plg_t *polygones, int npolygones);
extern int plg_cartesian(projPJ ref, vector<plg_t> & polygones);

extern frame_t plg_spherical_minmax(const plg_t *polygones, int npolygones);
extern frame_t plg_spherical_minmax(const vector<plg_t> &polygons);

extern int CheckPolar(vector<plg_t> polygons, projPJ &proj);
extern int CheckPolar(plg_t & p, projPJ &proj);

extern frame_t plg_cartesian_minmax(const plg_t & polygone);
extern frame_t plg_cartesian_minmax(const plg_t *polygones, int npolygones);
extern frame_t plg_cartesian_minmax(const vector<plg_t> & polygons);

extern  double plg_distance(const plg_t & plg_a, int vertex_a, const plg_t & plg_b, int vertex_b, int mode);
extern  double plg_distance(plg_t plg_a, int vertex_a, int vertex_b, int mode);


extern frame_t plg_recale(plg_t *polygones, int npolygones);
extern frame_t plg_recale(vector<plg_t> & polygons, double center=NAN);

extern  int plg_recale(vector<plg_t> & polygons, frame_t frame, int mode);
extern  int plg_recale(plg_t *polygones, int npolygones, frame_t frame, int mode);

extern void plg_degree_recale(vector<plg_t> *polygons, double center=0.);

extern int plg_flip(plg_t *polygones, int target);
extern int plg_flip(plg_t polygon);

extern int plg_format_from_name(const string &format_name);
extern string plg_format_extension(int format);
extern int plg_print_formats();

extern int plg_print(const plg_t & p);

extern int plg_find_format(const string &filename);

extern int plg_load(const string &filename,int format,plg_t **polygones,int *npolygones);
extern int plg_load(const string &filename,plg_t **polygones,int *npolygones);
extern int plg_load(const string &filename, int format, vector<plg_t> & polygons);
extern int plg_load(const string &filename, vector<plg_t> & polygons);

extern int plg_save(const char *filename, int format, const plg_t *polygones, int npolygones);
extern int plg_save(const char *filename, const plg_t & polygone);
extern int plg_save(const char *filename, int format, const vector<plg_t> & polygones);

extern int plg_identify(int id, plg_t *polygones, int npolygones);

//

extern int plg_colocalize(plg_t & p, int m, plg_t & q, int n, int option);

extern void plg_copy_point(plg_t *dest, int destPoint, const plg_t & src, int srcPoint);
extern void plg_copy_point(plg_t *polygone, int dest, int src);
extern int plg_move_point(plg_t & polygon, int position, const plg_point_t & point);
extern int plg_move_point(plg_t & dest, int destPoint, double x, double y, double t=NAN, double p=NAN);

extern int plg_insert_point(plg_t & polygon, int position, const plg_point_t & point);
extern int plg_insert_point(plg_t *polygones, int target, int position, const plg_point_t & point);

extern int plg_delete_point(plg_t *polygone, int point);

extern int *plg_increase_entries (plg_t **polygones, int *npolygones, int additional);
extern int *plg_increase_reuse_entries (plg_t **polygones, int *npolygones, int additional);

extern int plg_add(vector<plg_t> *array,const plg_t & plg);
extern int plg_add(vector<plg_t> *array,double x, double y, double t=NAN, double p=NAN);
//

extern int plg_merge (vector<plg_t> & p, vector<plg_t> q);

extern int plg_concat (plg_t & pp, plg_t & qq, int mode=0);
extern int plg_concat (int target[2], plg_t *polygones, int npolygones, int mode=0);
extern int plg_concat (int target[2], vector<plg_t> & polygones, int mode);
extern int plg_concat (int p1, int p2, vector<plg_t> & polygones, int mode);
extern int plg_concat (plg_target_t target[2], plg_t *polygones, int npolygones);

extern int *plg_cut_polygone (plg_t **polygones, int *npolygones, int target, int point);
extern int plg_cut_polygone(vector<plg_t> & polygones, int target, int point, int force=0);

extern int plg_swap_entries (plg_t *polygones, int m, int n);
extern int plg_swap_entries (vector<plg_t> polygones, int m, int n);

//

extern int plg_load_scan(const char *filename, plg_t **polygones, int *npolygones);
extern int plg_write_scan(const char *filename, const plg_t **polygones, int *npolygones);
extern int plg_write_pocformula(const char *filename, const plg_t **polygones, int *npolygones);

extern int plg_inquire_shp(const char *filename, int *np, int *nblocks, int *shapetype);
extern int plg_read_shp(const char *filename, plg_t *polygones, int np, int mode);
extern int plg_write_shp(const char *filename, const plg_t *polygones, int np);

extern int plg_inquire_netcdf(const char *filename, int *np);
extern int plg_read_netcdf   (const char *filename, plg_t *polygones, int np);
extern int plg_write_netcdf  (const char *filename, const plg_t *polygones, int np);

extern int plg_inquire_gshhs (const char *filename, int *ns, int *np, int *fmt);
extern int plg_read_gshhs    (const char *filename, plg_t *polygones, int fmt);

extern int plg_inquire_cst (const char *filename, int *ns, int *np);
extern int plg_read_cst    (const char *filename, plg_t *polygones);
extern int plg_write_cst   (const char *filename, const plg_t *polygones, int np);

extern int plg_inquire_xiso (const char *filename, int *ns, int *np);
extern int plg_read_xiso (const char *filename, plg_t *polygones);

extern plg_array_t   plg_readneigh (const string &polygonsFileName) throw ();
extern vector<plg_t> plg_readneigh (const string &polygonsFileName, bool);

extern int plg_inquire_boundaries (const char *filename, int *ns, int *np);
extern int plg_read_boundaries (const char *filename, plg_t *polygones);
extern int plg_write_boundaries (const char *filename, const plg_t *polygones, int npolygones);

extern int plg_read_descriptor (const char *filename, plg_t* & boundaries, int & nboundaries, bool check, bool debug);
extern int plg_read_descriptor (const char *filename, vector<plg_t> & boundaries, bool check, bool debug);

extern int plg_write_descriptor (string rootname, vector<plg_t> & boundaries, bool debug);
extern vector<plg_t> plg_split(const vector<plg_t> & p, bool debug);

/*------------------------------------------------------------------------------
  re-sampling */
extern plg_t plg_resample(const plg_t & input, int no);
extern plg_t plg_resample(const plg_t & target, double radius, int equalize, bool debug=false);
extern vector<plg_t> plg_resample(const vector<plg_t> & target, double radius, int equalize, bool debug);
extern plg_t plg_resample(const plg_t & target, const vector<double> & prescribed, int equalize, bool debug);

extern int plg_decimate(plg_t & polygon, double threshold, bool debug);
extern int plg_decimate(vector<plg_t> & polygons, double threshold, bool debug);
extern int plg_decimate(plg_t *polygone,int decimation, bool debug);

extern plg_t plg_sample_shorelines(double t1, double p1, double t2, double p2, plg_t *shorelines, int nshorelines, double radius, int equalize);

extern plg_t plg_subdivide(plg_t target, double radius, int equalize, int style, bool debug);
extern vector<plg_t> plg_subdivide(const vector<plg_t> & target, double radius, int equalize, int style, bool debug);

/*------------------------------------------------------------------------------
  "shorelines" sampling */
extern plg_t plg_resample_rigid_obsolete(plg_t target, projPJ projection, double radius,  int equalize, int verbose=1);
extern plg_t plg_resample_rigid(plg_t target, projPJ projection, double radius,  int equalize, int verbose=1);
extern plg_t plg_resample_rigid_standard(plg_t target, projPJ projection, double radius, char mask, int equalize, int verbose);
extern plg_t plg_resample_rigid_variable(const plg_t & target, projPJ projection,const SGfield_t<float> & resolution, double radius, char mask, int equalize, int verbose);

/*------------------------------------------------------------------------------
  "open limits" sampling using node density field*/
extern plg_t plg_sample(plg_t target, grid_t & grid, float *density, float mask, criteria_t & criteria,  int equalize, bool debug);

/*------------------------------------------------------------------------------
  academic shapes */

extern int plg_cartesian_ellipse(double x0, double y0, double a, double b, range_t<double> angles, int npt, plg_t & plg);
extern int plg_create_ellipse(plg_array_t *plg, double a, double b, double t0, double p0);

extern int plg_cartesian_circle(double x0, double y0, double radius, range_t<double> angles, int npt, plg_t & plg);
extern int plg_cartesian_circle(double t0, double p0, projPJ proj, double radius, range_t<double> angles, int npt, plg_t & plg);

extern int plg_circles(const plg_t & p, projPJ proj, double factor, vector<plg_t> & q, vector<double> & sizes);
extern int plg_circles(const vector<plg_t> & p, projPJ proj, double factor, vector<plg_t> & q, vector<double> & sizes);

/*------------------------------------------------------------------------------
   */
extern int plg_compact_entries (plg_t **polygones, int *npolygones);
extern int plg_compact_entries (vector<plg_t> & polygons);

extern int plg_add_entry(vector<plg_t> & polygons, plg_t p);
extern int plg_add_entries(vector<plg_t> & polygons, plg_t *p, int n);
extern int plg_add_entries(vector<plg_t> & polygons, vector<plg_t> &p);

extern vector<plg_t> plg_duplicate(const vector<plg_t> & p0);
extern int plg_duplicate(const vector<plg_t> & p, vector<plg_t> & q);

extern int plg_destroy_entries (vector<plg_t> & polygons);

/*------------------------------------------------------------------------------
   */
extern vector<plg_t> plg_FindCountours(const grid_t & zgrid, float*array, float mask, const vector<float> & values, mesh_t & mesh, float* & UGarray, bool show_mesh, bool debug, int verbose=2);
extern vector<plg_t> plg_FindCountours(const grid_t & zgrid, float*array, float mask, float value, mesh_t & mesh, float* & UGarray, bool show_mesh, bool debug, int verbose=2);

extern vector<plg_t> plg_smooth(const string & rootname, const vector<plg_t> & polygons, plg_filter_t plg_filter, vector<float> & isovalues, float target, point2D_t point, int mode, bool debug);
extern int plg_CreateLimits(const string & ShorelinesFile, const plg_t & polygon, string rootname, plg_filter_t & plg_filter, point2D_t point, double isocontour, bool smooth, bool removePools, bool debug);

extern plg_t plg_GridLimits(const grid_t & grid, int increment=1);
extern vector<plg_t>  plg_GridLimits(grid_t grid, float *mask, bool debug, int verbose);

extern vector<plg_t> plg_array2vector(plg_array_t polygons);
extern vector<plg_t> plg_array2vector(plg_t *polygons, int npolygons);

extern plg_array_t   plg_vector2array(vector<plg_t> polygons);
extern int plg_vector2array(vector<plg_t> polygons, plg_t* & p, int & np);

extern void plg_smoothlines(plg_t & p);
extern void plg_smoothlines(vector<plg_t> & p);

extern vector<plg_t> plg_extract(const string & ShorelinesFile,   const vector<plg_t> & selection, const string & rootname, int coordinates, int target, bool debug);
extern vector<plg_t> plg_extract(vector<plg_t> & shorelines_base, const vector<plg_t> & selection, const string & rootname, int coordinates, int target, bool debug);

extern int plg_RemoveInlandLimits(vector<plg_t> & p, double x, double y, projPJ projection, int check, bool debug);
extern int plg_RemoveImbrications(vector<plg_t> & p, double x, double y, int check);


extern int plg_export2nei(vector<plg_t> & p, mesh_t & mesh, bool strict, bool do_edges, bool debug);
extern  vector<plg_t> Nei2Polygons (mesh_t & mesh, bool check_islands, int verbose, bool debug) ;

extern int plg_reduce_selection(vector<plg_t> & polygons, const vector<plg_t> & selection, int mode, bool debug);
extern int plg_CutAtIntersections(vector<plg_t> & p, vector<plg_t> & q, int mode);

extern int plg_CutAtIntersectionsMultiples(const vector<plg_t> & p, vector<plg_t> & q, vector<plg_t> & islands, vector<plg_t> & external, vector<plg_t> & targeted, const vector<plg_t> & limits, bool debug);

extern int plg_SplitPolygons(vector<plg_t> & polygons, const vector<plg_t> & selection, vector<plg_t> & outside, vector<plg_t> & inside, vector<plg_t> & intersecting, bool strict, int coordinates);

extern int plg_QuickFrameReduction(vector<plg_t> & polygons, frame_t frame, int coordinates);
extern int plg_PolygonReduction(vector<plg_t> & polygons, const vector<plg_t> & selection, bool stric, int coordinates, int verbose, bool debug);

extern int plg_DeleteSelection(vector<plg_t> & polygons, const vector<plg_t> & selection, bool strict, int coordinates);

extern vector<plg_t> plg_extract_out(vector<plg_t> & shorelines_base, const vector<plg_t> & selection, string rootname, int mode, bool debug);


extern int plg_SetZone(const char *zone, frame_t & frame, string & rootname, point2D_t & point);
extern int plg_DecodeFrame(const string input, frame_t & frame);
extern int plg_DecodePosition(const string & input, point2D_t & point);


extern vector<plg_t>  plg_dilatation_cartesian(vector<plg_t> & p, double factor, bool debug);

extern int plg_SortBySize(vector<plg_t> & p, int verbose);
extern int plg_CheckClosure(vector<plg_t> & shorelines_base, bool repair, int verbose, bool debug);

extern int plg_CheckDuplicated(vector<plg_t> & polygons, bool repair, string rootname, bool verbose);
extern int plg_CheckDuplicated(plg_t *polygons, int npolygons, bool repair, string rootname, bool verbose);


extern projPJ plg_DefineProjection(plg_t *polygones, int npolygones);
extern projPJ plg_DefineProjection(vector<plg_t> & polygons, char *parmaters);
 
extern int plg_join (int target[2], plg_t *polygones, int npolygones, double join_radius=500.);
extern int plg_join (int target[2], vector<plg_t> polygones, double join_radius=500.);
extern int plg_join (plg_t & p, plg_t & q, double join_radius=500.);

extern float plg_resolution(plg_t & target, int k, grid_t & grid, float *density, float mask, criteria_t & criteria,  int equalize);

//extern plg_t plg_sample(plg_t target, grid_t grid, float *density, float mask, criteria_t criteria,  int equalize);

#endif  // #ifndef POLYGONES_H
