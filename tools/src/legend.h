
#if LEGEND_H == 0
#define LEGEND_H 1

#include "functions.h"

class point_t {
private :
public :
  int N;
  int nz;
  double t,p;
  float *z;
  point_t() {
    N=-1;
    t=p=1.e+35;
    z=0;
    nz=0;
    }
  point_t(int id, double x, double y, int ncols) {
    N=id;
    t=x;
    p=y;
    nz=ncols;
    z=new float[nz];
    }
  void destroy() {
    deletep(&z);
    N=-1;
    t=p=1.e+35;
    }
};

class point3D_t {
private :
public :
  int N;
  float t,p;
  float **z;
  float *depths;
  
  point3D_t() {
    N=-1;
    t=p=1.e+35;
    z=0;
    depths=0;
  }
};

class legend_t {
private :
public :
   int    ID, Type;
   double T,P,X,Y;
   char *ptr;
   
   legend_t() {
     ID=-1;
     Type=-1;
     T=P=X=Y=0.;
     ptr=0;
     }
     
//   void destroy() {
//     switch(this->Type) {
//       case LGD_TEXT_SYMBOL:
//         break;
//       case LGD_GRAPH:
//         legend02_t *tmp=(legend02_t *) this->ptr;
//         tmp.destroy();
//         break;
//       case LGD_PROFILE:
//         break;
//       }
//     }
  };

class legend00_t {
private :
public :
  int   ID, Type;
  float T,P,X,Y;
  int   I,J;
  int   coordinate_system;
  char  *text;
  int   fontID, FontType, FontSize;
  int   angle,  Hjustify, Vjustify;
  int   TextFG, TextBG;
  bool framed;
};

class legend01_t {
private :
public :
  int  ID, Type;
  float T,P,X,Y;
  int  I,J;
  char  Text[256];
  int  Font, FontType, FontSize;
  int  Orientation, Hjustify, Vjustify;
  int  TextFG, TextBG;
  bool TextFramed;
  int  Symbol;
  float SymbolSize;
  int  PenColour, FillColour;
  int  Xshift,Yshift;
  };

class legend02_t {
private :
public :
  int  ID,np,nz;
  int  currentz;
  point_t *points;
  int  pen0,pen1;
  int  framed;
  legend01_t *text;
  
  legend02_t() {
    ID=np=nz=-1;
    currentz=-1;
    points=0;
    text=0;
    }
    
  legend02_t(int npoints, int ncols) {
    np=npoints;
    nz=ncols;
    currentz=0;
    points=new point_t[np];
    for (int n=0; n < np;n++) {
      points[n].N=n;
      points[n].z=new float[nz];
      }
    text=0;
    }
    
  void destroy() {
    points->destroy();
    deletep(&text);
    ID=np=nz=-1;
    currentz=-1;
    }
  };

class legend03_t {
private :
public :
  int  ID,np,ndepths,ncolumns;
  int  currentz;
  point3D_t *points;
  int  pen0,pen1;
  int  framed;
  legend01_t *text;
  
  legend03_t() {
    ID=np=ndepths=ncolumns=-1;
    currentz=-1;
    points=0;
    text=0;
    }
};

/*#######################################################################
 ---- PARAMETERS -------------------------------------------------------
 #######################################################################*/


/*#######################################################################
 ---- VARIABLES --------------------------------------------------------
 #######################################################################*/

#ifdef  LEGEND_MAIN
#define PREFIX
#else
#define PREFIX extern
#endif

PREFIX int 	lgd_n,lgd_modified;
PREFIX int 	lgd_nmax;
PREFIX int 	lgd_lastedited;
PREFIX int 	lgd_defaulttype,lgd_graphformat;

PREFIX legend_t *legendstab;

#undef PREFIX

/*#######################################################################
 ---- PROTOTYPES --------------------------------------------------------
 #######################################################################*/

/* legend01C.c */
extern void lgd_free(int *);

extern void lgd_init(int *);

extern void lgd_alloc(int,int *);

extern void lgd_realloc(int,int *);

extern void lgd_clearall(int *);

extern void lgd_tocartesian(legend_t *, int *, int *);

extern int lgd_save(const char *out, legend_t *legend, int nlegend, char *fmt, char *comments);

extern void lgd_defaultlegend(legend_t *,int *);

extern void lgd_load(const char *in,  int *rstatus);

extern legend_t *lgd_import(double *x, double *y, float  **z, int np, int ncols);
extern legend_t *lgd_import(double *x, double *y, double **z, int np, int ncols);
extern legend_t *lgd_import(const vector<point_t> & points);

#endif
