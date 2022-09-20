


#ifndef ACADEMIC_H
#define ACADEMIC_H

#define TOPO_UNDEFINED    -1
#define TOPO_FIXED         0
#define TOPO_TRANSITION    1
#define TOPO_RESIZABLE     2
#define TOPO_SECTION_FILE  3

#define TOPO_LINEAR       0
#define TOPO_QUADRATIC    1
#define TOPO_EXPLICIT     2

#include "xyz.h"
#include "map.h"

class profile_t {
private :
public :
  int properties, type;
  double *x, *y, *d, *z, mask;
  int n;
  
  profile_t() {
    x=y=d=z=0;
    mask=NAN;
    n=-1;
    }
    
  void rescale(double s) {
    for(int k=0;k<n;k++) {
      if(z[k]!=mask) z[k]*=s;
      }
    }
    
  void flip() {
    for(int k=0;k<n/2;k++) {
      double tmp=z[k];
      z[k]=z[n-1-k];
      z[n-1-k]=tmp;
      }
    }
    
  double set_abscisse() {
    double size;
    deletep(&d);
    d=new double[n];
    d[0]=0;
    for(int k=1;k<n;k++) {
      d[k]=d[k-1]+geo_distance(x[k-1],y[k-1],x[k],y[k]);
      }
    size=d[n-1];
    return(size);
    }
    
  double interpolate(double dd) {
    double value;
    int status;
    status=map_interpolate1D(z, d, mask, n, dd, &value);
    return(value);
    }
};

class atom_shape_t {
private :
public :
  int properties, type;
  double length;
  double limits[2];
  double slope,h[10];
  double curving;
  string section_file;
  profile_t p;
 
  atom_shape_t() {
    properties=TOPO_UNDEFINED;
    length=slope=0;
    limits[0]=limits[1]=-1.;
    }

  void set_linear(int t, double l, double a, double h0, double h1) {
    properties=t;
    type=TOPO_LINEAR;
    if(isnan(l))  l=(h1-h0)/a;
    if(isnan(a))  a=(h1-h0)/l;
    if(isnan(h0)) h0=h1-a*l;
    if(isnan(h1)) h1=h0+a*l;
    length=l;
    h[0]=h0;
    h[1]=h1;
    slope=a;
    }

  void set_parabolic(int t, double l, double h0, double h1, double h2) {
    properties=t;
    type=TOPO_QUADRATIC;
    length=l;
    h[0]=h0;
    h[1]=h1;
    h[2]=h2;
    }

  void set_external(int t, string & s) {
    int status;
//     profile_t p;
    bool debug=false;
    properties=t;
    type=TOPO_EXPLICIT;
    section_file=s;
    status=xyz_loadraw (s.c_str(), (string) "YXZ", (char *) 0, p.x, p.y, p.z, &p.mask, p.n, debug);
    statistic_t ss=get_statistics(p.z, p.mask,p.n, 1);
    if(ss.mean<0.0) p.rescale(-1.0);
    p.flip();
    length=p.set_abscisse();
    }

  void reset() {
    properties=-1;
    }

  void destroy() {
    (*this).reset();
    }
    
  double actual_slope() {
    double s;
    s=(h[1]-h[0])/length;
    return(s);
    }

 int depth(double d, double & depth) {
    double x=0;
    int status=-1;
    double epsilon=1.e-07;
    
    d=d-limits[0];
    
    if(length>0) {
      x=d/length;
      if( (x>= 0.0-epsilon) && (x<=1.0+epsilon)) {
        switch (type) {
          case TOPO_LINEAR:
            depth=h[0]+slope*d;
            break;
          case TOPO_QUADRATIC:
            depth=2.0*((1.0-x)*(0.5-x)*h[0]+2.0*x*(1.0-x)*h[1]-x*(0.5-x)*h[2]);
            break;
          case TOPO_EXPLICIT:
            depth=p.interpolate(d);
            break;
          }
        status=0;
        }
      else {
        status=-1;
        }
      }
    else {
      status=-1;
      }

    return(status);
    }

};

class shape_t {
private :
public :
  int proportional,type;
  vector<atom_shape_t> transects;
  
  shape_t() {
    type=-1;
    proportional=0;
    }

  void reset() {
    type=-1;
    proportional=0;
    }

  void destroy() {
    for(int k=0;k<transects.size();k++) {
      transects[k].destroy();
      } 
    transects.clear();
    this->reset();
    }
    
  int depth(double x, double & depth) {
    double L=0;
    int status=-1;
    
    for(int k=0;k<transects.size();k++) {
      status=transects[k].depth(x,depth);
      if(status==0) {
        printf("%s : %d %lf %lf\n",__func__,k,x,depth);
        break;
        } 
      } 
    
    return(status);
    }

  double length(double & resizable) {
    double d, L=0;
    double depth;
    int i,status;
        
    resizable=0;
    
    for(int k=0;k<transects.size();k++) {
      L+=transects[k].length;
      if(transects[k].properties==TOPO_RESIZABLE) resizable+=transects[k].length;
      }  
      
    return(L);
    }

  void resize(double scale) {
    double d, L=0;
    double depth;
    int i,status;

    for(int k=0;k<transects.size();k++) {
      if(transects[k].properties==TOPO_RESIZABLE) transects[k].length*=scale;
      }
    }
   
  int parse(double size) {
    double d, L=0, resizable, scale;
    double depth;
    int i,k,status;
        
/*------------------------------------------------------------------------------
    adjust lengths */
    L=this->length(resizable);
    
    scale=(size+resizable-L)/resizable;
    
    this->resize(scale);
    
    L=this->length(resizable);

/*------------------------------------------------------------------------------
    build sections limits */
    k=0;
    transects[k].limits[0]=0.0;
    transects[k].limits[1]=transects[k].length;
    
    for(int k=1;k<transects.size();k++) {
      transects[k].limits[0]=transects[k-1].limits[1];
      transects[k].limits[1]=transects[k].limits[0]+transects[k].length;
      }
    
/*------------------------------------------------------------------------------
    adjust depths in transition sections */
    for(int k=1;k<transects.size()-1;k++) {
      if(transects[k].properties==TOPO_TRANSITION) {
        double hh[2];
        transects[k].depth(transects[k].limits[0],hh[0]);
        transects[k].depth(transects[k].limits[1],hh[1]);
        transects[k-1].depth(transects[k].limits[0],transects[k].h[0]);
        transects[k+1].depth(transects[k].limits[1],transects[k].h[1]);
        transects[k].slope=(transects[k].h[1]-transects[k].h[0])/transects[k].length;
        }
      } 

    return(0);
    }
    
};

//   inline shape_t operator*(const double x,const shape_t & a){
//     shape_t b;
//     return b;
//     }
// 
//   inline shape_t operator+(const shape_t & a, const shape_t & b){
//     shape_t c;
//     return c;
//     }

class atom_river_t {
private :
  shape_t shape;
public :
  int properties, type;
  double length, resolution;
  double limits[2];
  double wavenumber, phase;
  size_t n;
  double width[2];
  double datum[2];
  shape_t shapes[2];
  
  void reset() {
    properties=TOPO_UNDEFINED;
    length=0;
    limits[0]=limits[1]=-1.;
    }

  atom_river_t() {
    this->reset();
    }

  void set(int t, double l, double r, double k, double G) {
    properties=t;
    length=l;
    resolution=r;
    wavenumber=k;
    phase=G;
    n=NINT(l/r)+1;
    resolution=l/((double) n-1.0);
    }
    
  void set3D(double w1, double w2, double d1, double d2) {
    datum[0]=d1;
    datum[1]=d2;
    width[0]=w1;
    width[1]=w2;
    }

  void getwidth(vector<double> & ww) {
    double x,w;
    for(int i=0;i<n;i++) {
      x=i*resolution/length;
      w=(1.0-x)*width[0]+x*width[1];
      }
    }

  double interpolate(int i, double buffer[2]) {
    double x,w;
    x=i*resolution/length;
    w=(1.0-x)*buffer[0]+x*buffer[1];
    return(w);
    }

  double interpolate_cosine(int i, double buffer[2]) {
    double x,w,wmin,wmax,delta;
    if(buffer[0]<buffer[1]) {
      wmin=buffer[0];
      delta=buffer[1]-buffer[0];
      x=i*resolution/length;
      }
    else {
      wmin=buffer[1];
      delta=buffer[0]-buffer[1];
      x=1.0-i*resolution/length;
      }
    x=2.0*x;
    w=wmin+delta*(1.0-cos(2*x*M_PI/2.0))/2.0;
    return(w);
    }

  int compute_shape(int i, double L, double l, double & depth) {
    int status;
    double x,w,z[2];
    x=i*resolution/length;
    shapes[0].parse(L);
    shapes[1].parse(L);
    status=shapes[0].depth(l, z[0]);
    if(status!=0) return(status);
    status=shapes[1].depth(l, z[1]);
    if(status!=0) return(status);
    depth=(1.0-x)*z[0]+x*z[1];
    if(fabs(depth)<0.01) {
      status=0;
      }
    return(status);
    }
    
    
  void destroy() {
    shapes[0].destroy();
    shapes[1].destroy();
    shape.destroy();
    this->reset();
    }
};

class river_t {
private :
public :
  int properties, type;
  vector<atom_river_t> parts;
  size_t n;
  
  river_t() {
    }

  void set() {
    }
    
  void parse() {
    n=0;
    for(int k=0;k<parts.size();k++) {
      n+=parts[k].n;
      }
    }
    
  void destroy() {
    for(int k=0;k<parts.size();k++) {
      parts[k].destroy();
      }
    parts.clear();
    }
    
};




class topo_section_t {
private :
public :
  int    type;
  double length, slope, h0;
  double limits[2];
  string section;

  topo_section_t() {
    type=TOPO_UNDEFINED;
    length=slope=h0=0;
    limits[0]=limits[1]=-1.;
    }

  void set(int t, double l, double a, double h) {
    type=t;
    length=l;
    slope=a;
    h0=h;
    }

  void reset() {
    type=-1;
    }

  void destroy() {
    (*this).reset();
    }
    
  int depth_linear(double d, double & depth) {
    double x=0;
    int status;
    
    d=d-limits[0];
    
    if(length>0) {
      x=d/length;
      if( (x>= 0.0) && (x<=1.0)) {
        depth=h0+slope*d;
        status=0;
        }
      else {
        status=-1;
        }
      }
    else {
      status=-1;
      }
    return(status);
    }
    
  int depth(double d, double & depth) {
    int status;
    switch (type) {
      case TOPO_UNDEFINED:
        break;
      case TOPO_FIXED:
      case TOPO_TRANSITION:
      case TOPO_RESIZABLE:
        status=depth_linear(d, depth);
        break;
      case TOPO_SECTION_FILE:
        break;
      default: 
        break;
      }    
    return(status);
    }

};

class topo_t {
private :
public :
  int proportional,type;
  topo_section_t ridge,deep,shelfbreak,shelf,shorebreak,shore,beach;
//   vector<topo_section_t> transects;
  
  topo_t() {
    type=-1;
    proportional=0;
    }

  void reset() {
    type=-1;
//     transects.clear();
    }

  void destroy() {
    (*this).reset();
    }
    
  int depth(double d, double & depth) {
    double L=0;
    int status;
    list<topo_section_t>::iterator it;
    list<topo_section_t> list;
    
    list.push_back(beach);
    list.push_back(shore);
    list.push_back(shorebreak);
    list.push_back(shelf);
    list.push_back(shelfbreak);
    list.push_back(deep);
    list.push_back(ridge);
    
    it=list.begin();
    do {
      status=(*it).depth(d,depth);
      it++;
      } while( (status==-1) && (it != list.end()) );
    type=-1;
    list.clear();
    
    return(status);
    }

  double length(double & resizable) {
    double d, L=0;
    double depth;
    int i,status;
    vector<topo_section_t>::iterator it,previous;
    
    vector<topo_section_t> list;
    list.push_back(beach);
    list.push_back(shore);
    list.push_back(shorebreak);
    list.push_back(shelf);
    list.push_back(shelfbreak);
    list.push_back(deep);
    list.push_back(ridge);
    
    resizable=0;
    
    it=list.begin();
    do {
      L+=(*it).length;
      if((*it).type==TOPO_RESIZABLE) resizable+=(*it).length;
      it++;
      } while( it != list.end() );
    
    i=0;
    beach=list[i]; i++;
    shore=list[i]; i++;
    shorebreak=list[i]; i++;
    shelf=list[i]; i++;
    shelfbreak=list[i]; i++;
    deep=list[i]; i++;
    ridge=list[i]; i++;
    
    list.clear();
    return(L);
    }

  void resize(double scale) {
    double d, L=0;
    double depth;
    int i,status;
    vector<topo_section_t>::iterator it,previous;
    
    vector<topo_section_t> list;
    list.push_back(beach);
    list.push_back(shore);
    list.push_back(shorebreak);
    list.push_back(shelf);
    list.push_back(shelfbreak);
    list.push_back(deep);
    list.push_back(ridge);
        
    it=list.begin();
    do {
      if((*it).type==TOPO_RESIZABLE) (*it).length*=scale;
      it++;
      } while( it != list.end() );
      
    i=0;
    beach=list[i]; i++;
    shore=list[i]; i++;
    shorebreak=list[i]; i++;
    shelf=list[i]; i++;
    shelfbreak=list[i]; i++;
    deep=list[i]; i++;
    ridge=list[i]; i++;
    
    list.clear();
    }
   
  int parse(double size) {
    double d, L=0, resizable, scale;
    double depth;
    int i,status;
    vector<topo_section_t>::iterator it,previous,next;
    
    vector<topo_section_t> list;
    
/*------------------------------------------------------------------------------
    adjust lengths */
    L=this->length(resizable);
    
    scale=(size+resizable-L)/resizable;
    
    this->resize(scale);
    
    L=this->length(resizable);

    list.push_back(beach);
    list.push_back(shore);
    list.push_back(shorebreak);
    list.push_back(shelf);
    list.push_back(shelfbreak);
    list.push_back(deep);
    list.push_back(ridge);

/*------------------------------------------------------------------------------
    build sections limits */
    it=list.begin();
    status=(*it).limits[0]=0.0;
    status=(*it).limits[1]=(*it).length;
    do {
      previous=it;
      it++;
      status=(*it).limits[0]=(*previous).limits[1];
      status=(*it).limits[1]=(*previous).limits[1]+(*it).length;
      } while( it != list.end() );
    
/*------------------------------------------------------------------------------
    adjust depths in transition sections */
    it=list.begin();
    it++;
    do {
      if((*it).type==TOPO_TRANSITION) {
        previous=next=it;
        previous--;
        next++;
        double h[2];
        status=(*it).depth((*it).limits[0],h[0]);
        status=(*it).depth((*it).limits[1],h[1]);
        status=(*previous).depth((*it).limits[0],h[0]);
        status=(*next).depth    ((*it).limits[1],h[1]);
        (*it).h0=h[0];
        (*it).slope=(h[1]-h[0])/(*it).length;
        }
      it++;
      } while( it != list.end() );
    

    i=0;
    beach=list[i]; i++;
    shore=list[i]; i++;
    shorebreak=list[i]; i++;
    shelf=list[i]; i++;
    shelfbreak=list[i]; i++;
    deep=list[i]; i++;
    ridge=list[i]; i++;
    
    list.clear();
    return(0);
    }
}; 

extern int fe_spherical(mesh_t & mesh, projPJ projection);
extern int fe_create_Estuary(bool debug);


#endif
