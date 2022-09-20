
/*******************************************************************************

  T-UGO tools, 2006-2016

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief finite element clases
*/
/*----------------------------------------------------------------------------*/

#ifndef FE_CLASSES_H
#define FE_CLASSES_H

#include "functions.h" //safely includes omp.h
#include "discretisation.h"

class triangle_t {
private :
public :
  int vertex[3];                 /* vertex list                      */
  int edges[3];                  /* edges list                       */
  int ngh[3],nngh;               /* element neighbour list           */
  int opposite[3];               /* cross edge neighbour             */
  int *surrounding,nsurrounding; /* element surrounding list         */
  double DP[3], DQ[3];           /* Lagrange derivative coefficients */
  double nx[3], ny[3], l[3];     /* element side orientation/lengths */
  double dxdt, dydt;             /* element derivation coefficients  */
  double dxdp, dydp;             /* element derivation coefficients  */
  double dtdx, dtdy;             /* element derivation coefficients  */
  double dpdx, dpdy;             /* element derivation coefficients  */
  double jacobien;               /*   */
  double t_base,p_base;          /*   */
  double Area, TrueArea;         /* pseudo and true area             */
  double cosinus;                /*   */
  int child,ancestor,origin;    /* index in mesh hierarchy          */

  triangle_t() {
    nngh=-1;
    surrounding=0;
    nsurrounding=-1;
    child=ancestor=origin=-1;
    Area=TrueArea=NAN;
    }

  void destroy() {
    if(surrounding!=0) {
      delete[] surrounding;
      nsurrounding=0;
      };
    if(surrounding!=0) {
      delete[] surrounding;
      nsurrounding=0;
      };
    }
};

class vertex_t {
private :
public :
  double lon, lat, h;       /* longitude (degrees), latitude (degrees), mean depth */
  double c, s;              /* sine and cosine of latitude        */
  double *sigma,*zlevels;   /* sigma (0 at bottom, 1 at sea surface) and levels' immersion  */
  int *ngh, nngh;           /* neighbours (vertices) list         */
  int *edges,nedges;        /* connected edges list               */
  int *elmts,nelmts;        /* elements list                      */
  double mw;                /* mass weight (1/3 of elements area) */
  int code;                 /* boundary code for nodee            */
  int   child, ancestor;    /* index in mesh hierarchy            */
  int   origin;             /* index in mesh hierarchy            */

  void remove_neighbour(int N) {
    int i,j;
    for(i=0;i<nngh;i++) {
      if(ngh[i]==N) {
        for(j=i;j<nngh-1;j++)ngh[j]=ngh[j+1];
        ngh[nngh-1]=-1;
        nngh--;return;
        }
      }
    }
  
  void add_neighbour(int N) {
    int *vec_tab,i;
    nngh++;
    vec_tab=new int[nngh];
    for(i=0;i<nngh-1;i++)vec_tab[i]=ngh[i];
    vec_tab[nngh-1]=N;
    delete (ngh);
    ngh=vec_tab;
    }

  void null_value() {
    lon=9999.999;lat=9999.999; h=-1;
    sigma=NULL;ngh=NULL;nngh=-1;
    elmts=NULL;nelmts=-1;
    edges=NULL;nedges=-1;
    code=-1;ancestor=-1;
    }

  vertex_t() {
    lon=9999.999;lat=9999.999; h=-1;
    c=9999.999; s=9999.999;
    sigma=NULL;
    zlevels=NULL;
    ngh=  NULL;nngh=-1;
    elmts=NULL;nelmts=-1;
    edges=NULL;nedges=-1;
    code=-1;ancestor=-1;child=-1;
    }

  void allocate() {
    if(nngh>0)
      ngh=new int[nngh];
    if(nelmts>0)
      elmts=new int[nelmts];
    if(nedges>0)
      edges=new int[nedges];
    }
  
  void destroy() {
    deletep(&ngh);
    deletep(&edges);
    deletep(&elmts);
    }
};


class edge_t {
private :
public :
  int shared[2];        /* elements sharing this edge            */
  int eindex[2];              /* edge index in shared elements         */
  int vindex[2][2];     /* vertices index bounding the edge      */
  int nshared       ;   /* count of elements sharing this edge   */
  int extremity[2];     /* vertices bounding the edge            */
  int *ngh, nngh;       /* neighbours (edges) list               */
  int *vertex, nvertex;       /* neighbours (vertex) list              */
  double lon, lat;      /* edge centre coordinates               */
  double c, s;          /* edge centre cos/sin latitude          */
  double Tx, Ty, L;     /* tangent and length                    */
  double mw;            /* mass weight (1/3 of elements area)    */
  double *weight3D;     /* integration weight                    */
  int code;             /* boundary or interior edge             */
  int child,ancestor;   /* index in mesh hierarchy               */
  int origin;           /* index in mesh hierarchy               */

  edge_t() {
    lon=9999.999;lat=9999.999;
    c=9999.999;s=9999.999;
    ngh=   NULL;nngh=-1;
    vertex=NULL;nvertex=-1;
    weight3D=NULL;
    code=-1;ancestor=-1;child=-1;
    vindex[0][0]=vindex[0][1]=vindex[1][0]=vindex[1][1]=-1;
    eindex[0]=eindex[1]=-1;
    shared[0]=shared[1]=-1;
    }

  void destroy() {
    deletep(&ngh);
    deletep(&vertex);
    }
  
  int opposite_vertex(int n) {
    int nn=( (n==extremity[0]) ? extremity[1] : extremity[0]);
    return(nn);
    }
//   int opposite(int m) {
//     }

};

class face_t {
private :
public :
  int shared[2];      /* elements sharing this edge            */
  int nshared;        /* count of elements sharing this edge   */
  int extremity[2];   /* vertices bounding the edge            */
//   double Tx, Ty, L;   /* tangent and length                    */
//   double mw;          /* mass weight (1/3 of elements area)    */
//   int code;           /* boundary or interior edge             */
//   int ancestor;       /* index in parent mesh                  */
  face_t() {
    nshared=0;
    }

};

typedef struct {
  triangle_t *triangle;    /*                                  */
  double h[2][3];         /*                                  */
  double dxdt,dydt;       /* element derivation coefficients  */
  double dxdp,dydp;       /* element derivation coefficients  */
  double *weight;       /* integration weight               */
  int ancestor;           /* index in parent mesh             */
} prism_t;

class limit_t {
private :
public :
  int *vertex, nvertex; /* vertex index list                  */
  int *edges,  nedges;  /* edge index list                    */
  int code, *flags;     /* boundary code and type             */
  double length;        /* total limit length                 */
  double x,y;           /* hole point                         */

  void destroy() {
    if(vertex!=0) {
      delete[] vertex;
      vertex=0;
      };
    if(edges!=0) {
      delete[] edges;
      edges=0;
      };
    }

  limit_t() {
    vertex=0 ; nvertex=-1;
    edges=0;   nedges=-1;
    code=-1;   flags=0;
    length=-1;
    x=0; y=0;
    }
};

typedef struct {
  int *vertex;     /* vertex index list                  */
} neighbours_t;

#include "quadrangle.h"

class mesh_t {
private :
public :
  vertex_t     *vertices;           /* vertices set             */
  triangle_t   *triangles;          /* triangle element set     */
  quadrangle_t *quadrangles;        /* quadrangle element set   */
  edge_t       *edges;              /* edges list               */
  face_t       *faces;              /* faces list               */
  limit_t      *limits;             /* limits list              */
  int nvtxs, nedges, nfaces;        /*                          */
  int ntriangles, nquadrangles;     /*                          */
  int nelements;                    /*                          */
  int nlimits, nholes;              /*                          */
  int nnghm;                        /*                          */
  int **elt_nghl, *elt_nngh;        /* elements neigbour list   */
  int nlevels, nlayers;             /*                          */
  element_t   *elements;            /* for multi-element handling */
  hypermatrix_t matrix;             /* mass matrix              */
  int hbw;                          /* nodes half bandwidth     */
  int level, type, periodic;         /* mesh type*/
  int units, circular;              /* mesh type and degree     */
  mesh_t  *origin;                  /* mesh descriptor flag     */
  discretisation_t LGP0descriptor;  /* LGP0 mesh descriptor     */
  discretisation_t LGP1descriptor;  /* LGP1 mesh descriptor     */
  discretisation_t DGP1descriptor;  /* DGP1 mesh descriptor     */
  discretisation_t NCP1descriptor;  /* NCP1 mesh descriptor     */
  discretisation_t DNP1descriptor;  /* DNP1 mesh descriptor     */
  discretisation_t LGP2descriptor;  /* LGP2 mesh descriptor     */
  discretisation_t DGP2descriptor;  /* LGP2 mesh descriptor     */
  discretisation_t IPGdescriptor;   /* IPG mesh descriptor      */
  discretisation_t CQP0descriptor;  /* mesh descriptor flag     */
  discretisation_t CQP1descriptor;  /* mesh descriptor flag     */
  discretisation_t CQN1descriptor;  /* mesh descriptor flag     */
  
  int set_elements() {
#define FE_UNDEFINED    -1
#define FE_LINE          1
#define FE_TRIANGLE      2
#define FE_QUADRANGLE    3
#define FE_PENTACLE     20
    switch (nature()) {
      case FE_TRIANGLE:
        if(triangles==0) return(-1);
        nelements=ntriangles;
        elements=new element_t[nelements];
        for (size_t m=0;m<ntriangles;m++) {
          elements[m].vertex=triangles[m].vertex;
          elements[m].edges =triangles[m].edges;
          elements[m].size=3;
          elements[m].area=triangles[m].TrueArea;
//           elements[m].rank=triangles[m].rank;
//           elements[m].centralized=triangles[m].centralized;
//           elements[m].parent=triangles[m].parent;
          }
        break;
      case FE_QUADRANGLE:
        if(quadrangles==0) return(-1);
        nelements=nquadrangles;
        elements=new element_t[nelements];
        for (size_t m=0;m<nquadrangles;m++) {
          elements[m].vertex=quadrangles[m].vertex;
          elements[m].edges =quadrangles[m].edges;
          elements[m].size=4;
          elements[m].area=quadrangles[m].Area;
//           elements[m].rank=quadrangles[m].rank;
//           elements[m].centralized=quadrangles[m].centralized;
//           elements[m].parent=quadrangles[m].parent;
          }
        break;
      case FE_UNDEFINED:
        nelements=0;
        break;
      }
    return(0);
    }
  

  void destroy() {
    int n;
    if(vertices!=0) {
      for(n=0;n<nvtxs;n++) {
        vertices[n].destroy();
        }
      deletep(&vertices);
      };
    if(triangles!=0) {
      for(n=0;n<ntriangles;n++) {
        triangles[n].destroy();
        }
      deletep(&triangles);
      };
    if(quadrangles!=0) {
      for(n=0;n<nquadrangles;n++) {
        quadrangles[n].destroy();
        }
      deletep(&quadrangles);
      };
    if(edges!=0) {
      for(n=0;n<nedges;n++) {
        edges[n].destroy();
        }
      deletep(&edges);
      };
    deletep(&faces);
    if(limits!=0) {
      for(n=0;n<nlimits;n++) {
        limits[n].destroy();
        }
      deletep(&limits);
      };
    matrix.destroy();
    this->reset();
    LGP0descriptor.destroy();
    LGP1descriptor.destroy();
    DGP1descriptor.destroy();
    NCP1descriptor.destroy();
    LGP2descriptor.destroy();
    DGP2descriptor.destroy();
  }

private:
  void reset() {
    /* pointers */
    vertices=0;
    triangles=0;
    quadrangles=0;
    edges=0;
    faces=0;
    limits=0;
    origin=0;
    /* numbers */
    nvtxs=nedges=nfaces=0;
    ntriangles=nquadrangles=0;
    nelements=0;
    nlimits=nholes=0;
    nnghm=0;
    nlevels=nlayers=0;
    /* flags */
    hbw=-1;
    level=type=-1;
    periodic=0;
    units=circular=-1;
    }
public:
  
  mesh_t() {
    reset();
    }
  
  void allocate() {
    if(nvtxs>0)
      vertices=new vertex_t[nvtxs];
    if(ntriangles>0)
      triangles=new triangle_t[ntriangles];
    if(nquadrangles>0)
      quadrangles=new quadrangle_t[nquadrangles];
    if(nedges>0)
      edges=new edge_t[nedges];
    }
  
  int nature() const {
    if((ntriangles!=0) && (nquadrangles==0)) {
//      return(FE_TRIANGLE); /// HERE !!!
      return(2);
      }
    else if((ntriangles==0) && (nquadrangles!=0)) {
//      return(FE_QUADRANGLE); /// HERE !!!
      return(3);
      }
    return(-1);
    }
  
    mesh_t(string filename, int format);
    mesh_t& initialize(string filename, int format);
    mesh_t& initialize_descriptors(string FEdescription, int massmatrix = 1);
    mesh_t& initialize_descriptors(int FEdescription, int massmatrix = 1);

/*  int opposite_vertex(int m, int n) {
    int k[3];
    k[0]
    }*/
};

#define null_mesh mesh_t()

// class list_t {
// public :
//   int nx,ny;
//   vector<int> **elements;    ///< element indexes
//   omp_lock_t **locks;        ///< locks
//   double dx,dy;               ///<not used with mesh_t
//   
//   list_t(){
//     nx=0;
//     ny=0;
//     elements=NULL;
//     }
//   
//   void init(int gnx,int gny);
//   
//   ~list_t();
// };

class list_t {
public :
  int nx,ny;
  vector<int> **elements;    ///< element indexes
  omp_lock_t **locks;        ///< locks
  double dx,dy;               ///<not used with mesh_t
  
  list_t(){
    nx=0;
    ny=0;
    elements=NULL;
    locks=NULL;
    }
  
  void init(int gnx,int gny) {
  int k;
  
  nx=gnx;
  ny=gny;
  
  elements  =new vector<int> *[nx];
  locks     =new omp_lock_t  *[nx];
  exitIfNull(elements);
  
  for (k=0;k<nx;k++) {
    elements[k]=new vector<int> [ny];
    exitIfNull(elements[k]);
    locks[k]   =new omp_lock_t  [ny];
    exitIfNull(locks[k]);
    for(int l=0;l<ny;l++){
      elements[k][l].reserve(4);
      omp_init_lock(&locks[k][l]);
      }
    }
  }
  
  void destroy(void) {
    deletep2D(&locks,nx);
    deletep2D(&elements,nx);
    nx=0;
    ny=0;
    }
  
  ~list_t() {
    destroy();
    };
};

class numerics_t {
private :
public :
  int **nodes;
  int nnodes;
  int type;

  numerics_t() {
    nodes=0;
    nnodes=-1;
    type=-1;
    }
};

class location_t {
private :
public :
  char   *var;
  double lon, lat;
  char   *name;
  int    node;
  vertex_t  *vertex;
  
  location_t() {
    var=0;
    name=0;
  }
  void destroy() {
    deletep(&var);
    deletep(&name);
    }
};

typedef struct {
  edge_t edge;
  int node;
  int orientation;
} segment_t;

class section_t {
private :
public :
  location_t landmarks[2];
  segment_t  *segments;
  int nsegments;
  section_t() {
    segments=0;
    nsegments=-1;
  }
};

class sample_set_t {
private :
public :
  location_t *samples;
  int    n;
  sample_set_t() {
    samples=0;
    n=-1;
  }
  void destroy() {
    for(int k=0;k<n;k++) samples[k].destroy();
    deletep(&samples);
    }
};

class section_set_t {
private :
public :
  section_t *sections;
  int    n;
  section_set_t() {
    sections=0;
    n=-1;
  }
};

#endif
