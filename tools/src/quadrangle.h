

class quadrangle_t {
private :
public :
  int vertex[4];                    /* vertex list                            */
  int edges[4];                     /* edges list                             */
//   int *pNodes, *uNodes;             /* computational nodes list               */
//   int NpNodes, NuNodes;             /* computational nodes count              */
  int ngh[8],nngh;                  /* element neighbour list                 */
  int opposite[4];                  /* cross edge neighbour                   */
//  int *surrounding,nsurrounding;    /* element surrounding list               */
  double DP[4], DQ[4];              /* Lagrange derivative coefficients       */
  double nx[4], ny[4], l[4];        /* element side orientation/lengths       */
  double dxdt, dydt;                /* element derivation coefficients        */
  double dxdp, dydp;                /* element derivation coefficients        */
  double dtdx,dtdy;                 /* element derivation coefficients        */
  double dpdx,dpdy;                 /* element derivation coefficients        */
  double m0,m1,n0,n1,n2;            /* element localisation coefficients      */
  double J[4];                      /* element integration coefficients       */
  double jacobien;                  /*                                        */
  double Area, TrueArea;            /* pseudo and true area                   */
  double cosinus,tangente;          /*                                        */
  double t_base,p_base;             /*                                        */
  int child,parent,ancestro;        /* index in parent mesh                   */
  int origin;                       /* index in unpartitioned mesh            */
  int rank;                         /*                                        */

  quadrangle_t() {
//     pNodes=0;
//     uNodes=0;
//     NpNodes=-1;
//     NuNodes=-1;
    nngh=-1;
//     surrounding=0;
//     nsurrounding=-1;
    dxdt=dydt=dxdp=dydp=0;
    Area=TrueArea=-1;
    cosinus=tangente=0.;
    child=parent=ancestro=-1;
    origin=rank=-1;
    for(size_t k=0;k<4;k++) {
      vertex[k]=edges[k]=opposite[k]=-1;
      DP[k]=DQ[k]=0;
      nx[k]=ny[k]=l[k]=0;
      }
    for(size_t k=0;k<8;k++) {
      ngh[k]=-1;
      }
    }

  void destroy() {
//     if(pNodes!=0) {
//       delete[] pNodes;
//       pNodes=0;
//       };
//     if(uNodes!=0) {
//       delete[] uNodes;
//       uNodes=0;
//       };
//     if(surrounding!=0) {
//       delete[] surrounding;
//       nsurrounding=0;
//       };
    }

} ;

class element_t {
private :
public :
  int type;                      /* type of elements                 */
  int *vertex;                   /* vertex list                      */
  int *edges;                    /* edges list                       */
  int size;                      /* vertex/edge cardinal             */
  double area;                   /* element area                     */
  int centralized, parent, ancestor;       /* index in unpartitioned mesh      */
  int rank;                      /*                                  */

  element_t() {
    type=-1;
    vertex=0;
    edges=0;
    }

  void destroy() {
    }
};

class g_element_t {
private :
public :
  int *vertex;                   /* vertex list                      */
  int *edges;                    /* edges list                       */
  int nvertex,nedges;            /* edges list                       */
  double jacobien;               /*                                  */
  int child,parent,ancestro;     /* index in parent mesh             */
  int origin;                    /* index in unpartitioned mesh      */
  int rank;                      /*                                  */

  g_element_t() {
    vertex=0;
    edges=0;
    }

  void destroy() {
    }
} ;

class gmesh15_t {
private :
public :
  int  vertex;                   /* vertex list                      */

  gmesh15_t() {
    vertex=-1;
    }

  void destroy() {
    }
} ;

template <class Element> class gmesh_t {
private :
public :
  vertex_t  *vertices;             /* vertices set           */
  Element   *elements;             /* element set            */
  edge_t    *edges;                /* edges list             */
  face_t    *faces;                /* faces list             */
  limit_t   *limits;               /* limits list            */
  int **elt_nghl, *elt_nngh;       /* elements neigbour list */
  int nvtxs, nelts, nnghm, nedges, nlimits, nholes;   /*     */
  int nlevels, nlayers, nfaces;    /*                        */
  double *massM;                   /* mass matrix for LG P1  */
  int    *massP;                   /* mass pivot for LG P1   */
  hypermatrix_t matrix;            /* mass matrix            */
  int hbw;                         /* nodes half bandwidth   */
  int level, type, degree;         /* mesh type and degree   */
  int units, circular;             /* mesh type and degree   */
  int descriptor;                  /* mesh descriptor flag   */
  double lat_min,lat_max;          //ajout de thierry
  double lon_min,lon_max;          //ajout de thierry
  discretisation_t LGP0descriptor; /* LGP0 mesh descriptor   */
  discretisation_t LGP1descriptor; /* LGP1 mesh descriptor   */
  discretisation_t DGP1descriptor; /* DGP1 mesh descriptor   */
  discretisation_t NCP1descriptor; /* NCP1 mesh descriptor   */
  discretisation_t LGP2descriptor; /* LGP2 mesh descriptor   */

  void destroy() {
    int n;
    if(vertices!=0) {
      for(n=0;n<nvtxs;n++) {
        vertices[n].destroy();
        }
      delete[] vertices;
      vertices=0;
      };
    if(elements!=0) {
      for(n=0;n<nelts;n++) {
        elements[n].destroy();
        }
      delete[] elements;
      elements=0;
      };
    if(edges!=0) {
      for(n=0;n<nedges;n++) {
        edges[n].destroy();
        }
      delete[] edges;
      edges=0;
      };
    if(faces!=0) {
      delete[] faces;
      faces=0;
      };
    if(limits!=0) {
      for(n=0;n<nlimits;n++) {
        limits[n].destroy();
        }
      delete[] limits;
      limits=0;
      };
    matrix.destroy();
    (*this).reset();
    LGP2descriptor.destroy();
    DGP1descriptor.destroy();
    }

  void reset() {
    vertices=0;
    elements=0;
    edges=0;
    faces=0;
    limits=0;
    nvtxs=nelts=nnghm=nedges=nlimits=0;
    }

//   mesh_t<class Element>() {
//     vertices=0;
//     elements=0;
//     edges=0;
//     faces=0;
//     limits=0;
//     nvtxs=nelts=nnghm=nedges=nlimits=0;
//     }
};

