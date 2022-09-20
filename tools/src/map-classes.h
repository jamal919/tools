
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief structured grid classes declarations
*/
/*----------------------------------------------------------------------------*/

#if MAP_CLASSES_H == 0
#define MAP_CLASSES_H 1

#include <stdio.h>
#include <vector>
#include <string>
#include <iostream> /* for cout */

#include <proj_api.h> /* for projPJ */

#include "geo.h" /* for degree_recale */
#include "poc-time-classes.h" /* for date_t */
#include "fe-classes.h" /* for list_t */

#if TUGO
#include "functions.h" /* for range_t */
#else
#include "poc-array-stats.hpp" /* for range_t */
#endif


#define MODEH_UNSET -9
#define MODEH_REDUCED_GG -3
extern int ind1D(const int    *t,int n,int    time,int *m);

class grid_t
{
private:
public:
  int     nx, ny, nz, nt;
  double  dx, dy, dz;
  double  xmin, ymin, zmin;
  double  xmax, ymax, zmax;
  ///Number of dimensions of the x and y coordinates
  /**
  If 0, use min, max and d values.
  Can otherwise be 1 or 2.
  If -2: transposed.
  If -1: unstructured.
  */
  int modeH;
  ///Number of dimensions of the z coordinates
  /** Can only be 1 or 3 */
  int modeV;
  double  *x,*y,*z,*sigma,*bottom,zmask,*time;
  int *reduced_nx,*reduced_xi;/* for reduced gaussian grids */
  list_t *list;
  date_t  origine;
  signed char    *mask;
  double  dmask;
  int     circular,overlapped,connex;
  projPJ projection;
  char *proj4options;

  grid_t() {
    this->init();
    }

  void init() {
    nx=0; ny=0; nz=0; nt=0;
    dx=0; dy=0; dz=0;
    xmin=0; ymin=0; zmin=0;
    xmax=0; ymax=0; zmax=0;
    modeH=MODEH_UNSET;modeV=0;
    x=NULL;y=NULL;z=NULL;sigma=NULL;bottom=NULL;zmask=0;time=NULL;
    reduced_nx=NULL;reduced_xi=NULL;
    list=NULL;
    mask=NULL;
    dmask=0;
    circular=0;overlapped=0;connex=0;
    projection=0;
    proj4options=0;
    }
  
  void init(size_t nnx, size_t nny) {
    init();
    nx=nnx;
    ny=nny;
    x=new double[nx*ny];
    y=new double[nx*ny];
    }
  
  void xy(size_t i, size_t j, double & x, double & y) const {
  /// \sa map_grid_x(), map_grid_y()
    switch(modeH) {
      case 0:
        x=this->xmin+i*this->dx;
        y=this->ymin+j*this->dy;
        break;
      case 1:
        x=this->x[i];
        y=this->y[j];
        break;
      case 2:
      case -2:{
        size_t m=this->nx*j+i;
        x=this->x[m];
        y=this->y[m];
        }break;
      default:
        TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when modeH=%d\n",__func__,modeH);
      }
    }

  void set_xy(size_t i, size_t j, double x, double y){
    switch(modeH) {
      case 0:
        break;
      case 1:
        this->x[i]=x;
        this->y[j]=y;
        break;
      case 2:
      case -2:{
        size_t m=this->nx*j+i;
        this->x[m]=x;
        this->y[m]=y;
        }break;
      default:
        TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when modeH=%d\n",__func__,modeH);
      }
    }

  void ij(size_t m, int * i, int * j) const {
    if(modeH!=MODEH_REDUCED_GG){
      *i=m % this->nx;
      *j=m/this->nx;
      }
    else{
      int status,j_=-1;
      status=ind1D(reduced_xi,ny,m,&j_);
      if(reduced_xi[j_]>m)
        j_--;
      *i=m-reduced_xi[j_];
      *j=j_;
      }
    }

  void xy(size_t m, double & x, double & y) const {
    int i,j;
#if 0
    if(m>Hsize())
      TRAP_ERR_EXIT(ENOEXEC,"avoiding a crash");
#endif
    switch(modeH) {
      case 0:
        ij(m,&i,&j);
        x=this->xmin+i*this->dx;
        y=this->ymin+j*this->dy;
        break;
      case 1:
        ij(m,&i,&j);
        x=this->x[i];
        y=this->y[j];
        break;
      case 2:
      case -2:
        x=this->x[m];
        y=this->y[m];
        break;
      default:
        TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when modeH=%d\n",__func__,modeH);
      }
    }

  void set_xy(size_t m,const double x,const double y){
    int i,j;
    switch(modeH) {
      case 0:
        break;
      case 1:
        ij(m,&i,&j);
        this->x[i]=x;
        this->y[j]=y;
        break;
      case 2:
      case -2:
        this->x[m]=x;
        this->y[m]=y;
        break;
      default:
        TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when modeH=%d\n",__func__,modeH);
      }
    }

  double GetZ(size_t i, size_t j, size_t k) const {
    double z;
    switch(modeV) {
      case 0:
        z=this->zmin+k*this->dz;
        break;
      case 1:
        z=this->z[k];
        break;
      case 3:{
        size_t m=this->nx*(this->ny*k+j)+i;
        z=this->z[m];
        }break;
      default:
        TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when modeV=%d\n",__func__,modeV);
      }
    return(z);
    }

  double GetZ(size_t m, size_t k) const {
    switch(modeV) {
      case 0:
        return this->zmin+k*this->dz;
      case 1:
        return this->z[k];
      case 3:
        return this->z[this->ny*this->nx*k+m];
      default:
        TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when modeV=%d\n",__func__,modeV);
      }
    TRAP_ERR_EXIT(ENOEXEC,"coding error\n");
    return NAN;
    }

  void free() {
    deletep(&x);
    deletep(&y);
    deletep(&z);
    free_list();
    deletep(&sigma);
    deletep(&bottom);
    deletep(&time);
    deletep(&mask);
    nx=0; ny=0; nz=0; nt=0;
    deletep(&proj4options);
    deletep_void_void(&projection,pj_free);
    }

  void free_list() {
    if(list!=0) {
      list->destroy();
      deletep(&list);
      }
    }

  size_t Hindex(int i, int j) const {
    size_t n=(size_t) nx*(size_t) j+i;
    return(n);
    }

  size_t Hsize() const {
    size_t s=(size_t) nx*(size_t) ny;
    return(s);
    }

  size_t Vsize() {
    size_t s;
    switch(modeV) {
      case 0:
        s=0;
        break;
      case 1:
        s=this->nz;
        break;
      case 2:
        s=this->nz*this->ny*this->nx;
        break;
      default:
        s=-1;
        break;
        
      }
    return(s);
    }

  size_t xsize() const {
    size_t s;
    switch (modeH) {
      case 0:
        s=0;
        break;

      case 1:
        s=nx;
        break;

      case 2:
      case -2:
        s=nx*ny;
        break;

      default:
        TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when modeH=%d\n",__func__,modeH);
      }
    return(s);
    }

  size_t ysize() const {
    size_t s;
    switch (modeH) {
      case 0:
        s=0;
        break;

      case 1:
        s=ny;
        break;

      case 2:
      case -2:
        s=nx*ny;
        break;

      default:
        TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when modeH=%d\n",__func__,modeH);
      }
    return(s);
    }

//    void minmax() {
//     size_t s;
//     switch (modeH) {
//       case 0:
//         break;
// 
//       case 1:
//         range_t r=minmax(grid.x, grid.nx);
//         break;
// 
//       case 2:
//       case -2:
//         s=nx*ny;
//         break;
// 
//       default:
//         TRAP_ERR_EXIT(ENOEXEC,"programming error: %s called when modeH=%d\n",__func__,modeH);
//       }
//     return(0);
//     }

  size_t size() const {
    size_t s=(size_t) nx*(size_t) ny*(size_t) nz;
    return(s);
    }
  
  void duplicate(grid_t source) {
    *this=source;
    if(source.x!=0) {
      this->x=new double[source.xsize()];
      memmove(this->x,source.x, source.xsize()*sizeof(double));
      }
    if(source.y!=0) {
      this->y=new double[source.ysize()];
      memmove(this->y,source.y, source.ysize()*sizeof(double));
      }
    }

  void print(ostream &out=cout,const char *indent="") const {
    const int oldprecision=out.precision();
    out.precision(17);
    out <<indent<< "nx=" << nx << " ny=" << ny << " nz=" << nz << " nt=" << nt << endl;
    out <<indent<< "dx=" << dx << " dy=" << dy << " dz=" << dz << endl;
    out <<indent<< "xmin=" << xmin << " ymin=" << ymin << " zmin=" << zmin << endl;
    out <<indent<< "xmax=" << xmax << " ymax=" << ymax << " zmax=" << zmax << endl;
    out <<indent<< "modeH=" << modeH << " modeV=" << modeV << endl;
    out <<indent<< "zmask=" << zmask << " dmask=" << dmask << endl;
    out <<indent<< "circular=" << circular << " overlapped=" << overlapped << " connex=" << connex << endl;
    out.precision(oldprecision);
    out.flush();
    }

  void brief_print(ostream &out=cout,const char *indent="") const {
    out <<indent
      << xmin << "<=x<=" << xmax << " dx=" << dx << " nx=" << nx << " "
      << ymin << "<=y<=" << ymax << " dy=" << dy << " ny=" << ny << endl;
    out.flush();
    }

};

class Cgrid_t
{
private:
public:
  grid_t *z_grid, *f_grid, *u_grid, *v_grid;
  int *z_incidence, *f_incidence, *u_incidence, *v_incidence;
  Cgrid_t() {
    z_incidence = NULL;
    f_incidence = NULL;
    u_incidence = NULL;
    v_incidence = NULL;
    }
};

class resize_t
{
private:
public:
  vector<int> start, end, incr;
  vector<int> former;
  bool neutral;
  resize_t() {
    neutral=true;
    }
  void init(grid_t grid, int offset_i_low, int offset_i_high, int offset_j_low, int offset_j_high)
    {
    former.push_back(grid.nx);
    former.push_back(grid.ny);
    start.push_back(offset_i_low);
    start.push_back(offset_j_low);
    end.push_back(grid.nx-1-offset_i_high);
    end.push_back(grid.ny-1-offset_j_high);
    neutral=false;
    }
};

class gnames_t
{
private:
public:
  char *vlon, *vlat, *vtopo, *vmask;
  char *ha, *hg, *ua, *ug, *va, *vg;
  char *LSAa, *LSAg;
  gnames_t(){
    this->vlon=0;
    this->vlat=0; /* set but never used */
    this->vtopo=0;
    this->vmask=0;
    ha=hg = ua=ug = va=vg =0;
    LSAa=LSAg =0;
    };
  ~gnames_t(){
    deletep(&vlon);
    deletep(&vlat);
    deletep(&vtopo);
    deletep(&vmask);
    
    deletep(&ha);
    deletep(&hg);
    deletep(&ua);
    deletep(&ug);
    deletep(&va);
    deletep(&vg);
    
    deletep(&LSAa);
    deletep(&LSAg);
    };
};

class metagrid_t
{
private:
public:
  Cgrid_t Cgrid;
  gnames_t z_gnames, f_gnames, u_gnames, v_gnames;
  string format, gridfile, maskfile, topofile, target;
  string projection;
  float tag;
};


class frame_t {
public :
  double xmin,xmax,ymin,ymax;
  void init() {
    xmin=+INFINITY;
    xmax=-INFINITY;
    ymin=+INFINITY;
    ymax=-INFINITY;
    }
  frame_t() {
    init();
    }
  void init(double xxmin,double xxmax,double yymin,double yymax) {
    xmin=xxmin;
    xmax=xxmax;
    ymin=yymin;
    ymax=yymax;
    }
  frame_t(double xxmin,double xxmax,double yymin,double yymax) {
    init(xxmin,xxmax,yymin,yymax);
    }
  frame_t(const range_t<double> & x,const range_t<double> & y) {
    init(x.min,x.max,y.min,y.max);
    }
  frame_t(const grid_t & grid) {
    if(grid.modeH==MODEH_UNSET)
      init();
    else
      init(grid.xmin,grid.xmax,grid.ymin,grid.ymax);
    }
  int initialised() const {
    int status;
    if(isinf(xmin) and isinf(xmax) and isinf(ymin) and isinf(ymax)) {
      return(0);
      }
    status=(xmin<=xmax and ymin<=ymax);
    return status;
    }
  
  int inside(double x,double y) const {
    int status;
    status=((xmin<=x) && (x<=xmax) && (ymin<=y) && (y<=ymax));
    return(status);
    }
  
  double x_size() const {
    double size;
    if(xmin==+INFINITY or xmax==-INFINITY) size=0;
    else size=xmax-xmin;
    return(size);
    }
  
  double y_size() const {
    double size;
    if(ymin==+INFINITY or ymax==-INFINITY) size=0;
    else size=ymax-ymin;
    return(size);
    }
  
  double size() const {
    double size;
    size=this->x_size()*this->y_size();
    return(size);
    }
  
  double x_center() const {
    double size;
    size=0.5*(xmax+xmin);
    return(size);
    }
  
  double y_center() const {
    double size;
    size=0.5*(ymax+ymin);
    return(size);
    }
  
  void dilatation(double f) {
    double delta;
    delta=f*x_size()/2.0;
    xmin-=delta;
    xmax+=delta;
    delta=f*y_size()/2.0;
    ymin-=delta;
    ymax+=delta;
    }
  
  void dilatation(const frame_t f) {
    if(isfinite(xmin) and isfinite(xmax)){
      const double center=this->x_center();
      updatemin(&xmin,degree_recale(f.xmin,center));
      updatemax(&xmax,degree_recale(f.xmax,center));
      }
    else{
      updatemin(&xmin,f.xmin);
      updatemax(&xmax,f.xmax);
      }
    updatemin(&ymin,f.ymin);
    updatemax(&ymax,f.ymax);
    }
  
  frame_t intersection(const frame_t & f) const {
    frame_t intersection=*this;
    updatemax(&intersection.xmin,f.xmin);
    updatemin(&intersection.xmax,f.xmax);
    updatemax(&intersection.ymin,f.ymin);
    updatemin(&intersection.ymax,f.ymax);
    if(intersection.xmin>intersection.xmax) {
      intersection.xmin=+INFINITY;
      intersection.xmax=-INFINITY;
      }
    if(intersection.ymin>intersection.ymax) {
      intersection.ymin=+INFINITY;
      intersection.ymax=-INFINITY;
      }
    return(intersection);
    }
  bool intersect(const frame_t & f) const {
    if( f.inside(xmin,ymin) )
      return true;
    if( f.inside(xmax,ymax) )
      return true;
    if( f.inside(xmin,ymax) )
      return true;
    if( f.inside(xmax,ymin) )
      return true;
    if( inside(f.xmin,f.ymin) )
      return true;
    if( inside(f.xmax,f.ymax) )
      return true;
    return false;
    }
  void dilatation(double x,double y) {
    updatemin(&xmin,x);
    updatemax(&xmax,x);
    updatemin(&ymin,y);
    updatemax(&ymax,y);
    }
  void print(const char *s) {
    if(s!=0) printf("%s ",s);
    printf("xmin=%lf xmax=%lf ymin=%lf ymax=%lf\n",xmin,xmax,ymin,ymax);
    }
};

#endif
