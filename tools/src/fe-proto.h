/** \file
\brief prototypes of finite element functions
\date reviewed 1 Aug 2011
\author  Damien Allain

Link with :
-lfe -lrutin
*/

/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#ifndef FE_PROTO
#define FE_PROTO


#include "functions.h"
#include "polygones.h"

#include "legend.h"

#define ANSI_DECLARATORS
#ifndef __TRIANGLE
#include "triangle.h"
#define __TRIANGLE
#endif

#include "fe.def"

/* *-----------------------------------------------------------------------------
Basis functions */
extern void fe_LGP1base(double x, double y, double *beta);
extern void fe_LGP2base(double x, double y, double *beta);
extern void fe_LGP3base(double x, double y, double *beta);
extern void fe_LGP4base(double x, double y, double *beta);

extern void fe_QP0base(double x, double y, double *beta);
extern void fe_QP1base(double x, double y, double *beta);

extern int fe_LGPbase(const mesh_t & mesh,int discretisation,double x,double y,double *beta);

extern void fe_LGP1prime(double x, double y, double *beta_x, double *beta_y);

extern void fe_NCP1base(double x, double y, double *beta);

extern void fe_NCP1prime(double x, double y, double *beta_x, double *beta_y);
extern void fe_LGP2prime(double x, double y, double *beta_x,double *beta_y);

extern void fe_BL3Dbase (double x, double y, double z, double *beta);
extern void fe_BL3Dprime(double x, double y, double z, double *prime[3]);

extern int fe_gradient(mesh_t, discretisation_t, float *, int, float *,float *);
extern int fe_gradient(mesh_t, discretisation_t, double *, int, double *,double *);
extern int fe_gradient(mesh_t, discretisation_t, complex<float>  *, int, complex<float>  *, complex<float>  *);
extern int fe_gradient(mesh_t, discretisation_t, complex<double> *, int, complex<double> *, complex<double> *);

#include "fe-integration.h"

extern int fe_list(mesh_t *mesh) ;  
extern int fe_list(mesh_t *mesh, int ElementType);

extern int fe_e2n (mesh_t *mesh,int verbose=0);
extern int fe_element2neighbours (mesh_t mesh, discretisation_t descriptor);

extern int fe_bandwidth (mesh_t *mesh);
extern int fe_element_bw(mesh_t *mesh);

extern int fe_ReallocateVertices(mesh_t *, int);

extern int fe_geometry(mesh_t *mesh);

//avoid calling fe_initaffine() directly, this is most likely redundant
extern int fe_initaffine(mesh_t *mesh,int i, bool ignore_CW=true);
extern int fe_initaffine(mesh_t *mesh);

extern int fe_initaffine_spherical(const mesh_t & mesh, quadrangle_t *q, int m);

extern int fe_affine_directe(const triangle_t & e,double t,double p,double *ksi,double *eta,int mode);
extern int fe_affine_directe3D(prism_t prism, double t, double p, double h,double *x,double *y, double *z,int mode);

extern int fe_affine_inverse(triangle_t e,double *t, double *p,double ksi,double eta);

extern void fe_bilinear3D(double x,double y,double z,double *beta);

extern int fe_nearest_vertex(const mesh_t & mesh,double x,double y,int *except,int nexcept);
extern int fe_nearest_vertex(const mesh_t & mesh, double x, double y, int code, double & d);

extern int fe_nearest_node(const discretisation_t & descriptor,double x,double y,int *except,int nexcept);
  
extern int fe_nearest_boundaryvertex(const mesh_t & mesh,double x,double y,int *except,int nexcept);

extern int fe_whichnearest(const mesh_t & mesh,double t,double p,double *tnew,double *pnew);
extern int fe_testelement_v1(const mesh_t & mesh,int i,float t,float p,int mode);

extern int fe_testelement_v2(const mesh_t & mesh,const triangle_t & q,double t,double p,int mode=0);
extern int fe_testelement_v2(const mesh_t & mesh,const quadrangle_t & q,double t,double p,int mode=0);

extern int fe_testelement_v3(const mesh_t & mesh,const triangle_t & q,double t,double p);
extern int fe_testelement_v3(const mesh_t & mesh,const quadrangle_t & q,double t,double p);

extern int fe_whichelement_v2(const mesh_t & mesh,double t,double p);

extern int fe_whichelement(const mesh_t & mesh,double t,double p);
extern int fe_whichelement(const mesh_t & mesh,double t,double p, int hint);

extern int fe_whichelement_inlist(const mesh_t & mesh,const int *list,int nlisted,double t,double p);
extern int fe_whichelement_inlist(const mesh_t & mesh,const vector<int> & list,double t,double p);
extern int fe_whichelement_inlist(const mesh_t & mesh,const vector<int> & list,int nlisted,double t,double p);

extern void fe_AllocateList(const mesh_t & mesh,grid_t *grid,list_t *list);
extern void fe_CreateTriangleList(const grid_t & grid,const mesh_t & mesh,list_t *list);
extern void fe_UpdateList(const grid_t & grid,list_t *list,const range_t<double> &pR,range_t<double> tR, int i, int verbose=0);
extern void fe_Allocate_and_CreateList(const mesh_t & mesh,grid_t *grid, list_t *triangles, list_t *quadrangles=0);

extern int fe_InExtent(const mesh_t & mesh, const triangle_t   & element, double x, double y);
extern int fe_InExtent(const mesh_t & mesh, const quadrangle_t & element, double x, double y);

extern int fe_beta(mesh_t & mesh, double t, double p,int *elt,int *node, double *beta);
extern int fe_beta(mesh_t & mesh, double t, double p,int hint, int *elt,int *node, double *beta);

extern int fe_beta(mesh_t & mesh, const discretisation_t & descriptor, double t, double p, int *elements, int *nodes, double *beta);
extern int fe_beta(mesh_t & mesh, const discretisation_t & descriptor, double t, double p, int hint, int *elements, int *nodes, double *beta);

extern int fe_beta_inlist(const mesh_t & mesh,const int *list,int nlisted,double t,double p,int *elt,int *node,double *beta);
extern int fe_beta_inlist(const mesh_t & mesh,const vector<int> & list,double t,double p,int *elt,int *node,double *beta);

extern int fe_readmesh_QUODDY (const char *, const char *, const char *, mesh_t *, bool debug);
extern int fe_readmesh_TGL (const char *, const char *, mesh_t *);
extern int fe_readmesh_NGH (const char *,mesh_t *);

extern int fe_readmesh(const char *, int, mesh_t *, const char *proj4=0);
extern int fe_readmesh_TELEMAC_BINARY (const char *,mesh_t *, const char *);

extern int fe_ascii_saver1(const char *filename, mesh_t & mesh, float  *buffer, float  mask, int discretisation, int fmt, char **comments, int gArchiveExtent=0);
extern int fe_ascii_saver1(const char *filename, mesh_t & mesh, double *buffer, double mask, int discretisation, int fmt, char **comments, int gArchiveExtent=0);

extern int fe_ascii_saver2(const char *filename, mesh_t & mesh, double *bufx, double *bufy, double mask, int discretisation, int fmt, char **comments);
extern void fe_boundarycode(const char *filename, mesh_t & mesh, int level);

extern int fe_read_boundarycode (const char *, mesh_t & ,int);
extern int fe_write_boundarycode(const char *, const mesh_t &,int);
extern int fe_write_boundarycode(const char *belfile, vector<plg_t> & boundaries);

extern int fe_SaveLimits(const char *, const mesh_t &);

extern int fe_readnodes(const char *filename,int,mesh_t *mesh);
extern int fe_readnodes_TGL(const char *filename,mesh_t *mesh);

extern int fe_savenodes(const char *filename,int,mesh_t & mesh);
extern int fe_savenodes_TGD(const char *filename,mesh_t & mesh);
extern int fe_savepoly(const char *filename,mesh_t mesh);

extern int fe_savemesh_NGH(const char *filename,const mesh_t & mesh);
extern int fe_find_format(const string &filename);
extern int fe_savemesh(const char *, int, mesh_t &, const char *proj4=0);

extern int fe_updatemeshNC3D(const char *filename, mesh_t mesh, int option);

extern int fe_barycentre(mesh_t mesh,int i,double *t,double *p,int mode);

extern int fe_position(const mesh_t &mesh, const vertex_t &,double *t,double *p,int mode=0);
extern int fe_position(const mesh_t &mesh, const triangle_t &,double *t,double *p,int mode=0);
extern int fe_position(const mesh_t &mesh, const quadrangle_t & quadrangle, double *t,double *p,int mode=0);
extern int fe_position(const mesh_t &mesh,int m, double *t,double *p,int mode=0);
extern int fe_position(const mesh_t &mesh, const edge_t &,double *t,double *p,int mode=0);

extern void fe_integrale01(mesh_t mesh, float *buffer, double *sum);
extern void fe_integrale02(mesh_t mesh, float *buffer, double *sum, double latmax);
extern void fe_massmatrix01(mesh_t mesh, double **A, double*b, int ngauss);
extern void fe_massmatrix02(mesh_t mesh, double *A, int *pivot, int hbw);
extern void fe_massmatrix03(mesh_t mesh, double *A, int *pivot, int hbw, int ngauss);

extern const vector<int> *findInList(const list_t & list,const grid_t & grid,double dx,double dy,double t,double p,double *x,double *y,int verbose=0);

extern int fe_detectregular(const mesh_t &mesh,const grid_t &grid,int *buffer,float *fullness,int mode=0);
extern int *fe_detect(const mesh_t & mesh,const double *t,const double *p,int npositions);

extern int *fe_detect(const mesh_t & mesh,const triangle_t *,  const list_t & list,const grid_t & grid,const double *t,const double *p,int npositions);
extern int *fe_detect(const mesh_t & mesh,const quadrangle_t *,const list_t & list,const grid_t & grid,const double *t,const double *p,int npositions);

extern int fe_detect(const mesh_t & mesh,const triangle_t *,  const list_t & list,const grid_t & grid,double t,double p,int & last);
extern int fe_detect(const mesh_t & mesh,const quadrangle_t *,const list_t & list,const grid_t & grid,double t,double p,int & last);


extern int *fe_scan_elements(const mesh_t & mesh,const grid_t & grid,int mode=0,int option=0);
extern int *fe_scan_elements3d(mesh_t mesh,grid_t grid3d, int mode);

extern int fe_map(const mesh_t & mesh, float *buffer, int discretisation, const grid_t & grid, int *elts, float *c, float mask);
extern int fe_map(const mesh_t & mesh, double *buffer, int discretisation, const grid_t & grid, int *elts, double *c, double mask);
extern int fe_map(const mesh_t & mesh,complex<float> *buffer,int discretisation,const grid_t & grid,int *elts,complex<float> *c,complex<float> mask);
extern int fe_map(const mesh_t & mesh, complex<double> *buffer, int discretisation, grid_t grid, int *elts, complex<double> *c, complex<double> mask);

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int fe_map(const mesh_t & mesh,T *buffer,const grid_t & grid,int *elts,T *c,T mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=fe_map(mesh,buffer,LGP1,grid,elts,c,mask);
  
  return result;
}

extern int tide_atlasUG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, complex<double> *z,complex<double> cmask,const mesh_t &mesh,int discretisation,int verbose=1,const char *wave=0);
extern int tide_atlasUG2positions(const char *filename,const char *v1,const char *v2,const double *lon,const double *lat, int npositions, complex<float> *z,complex<float> cmask,const mesh_t &mesh,int discretisation,int verbose=1,const char *wave=0);

extern int fe_mapVdepth(mesh_t & mesh, grid_t grid,  short *c, short mask, int algorithm, int verbose);
extern int fe_mapVdepth(mesh_t & mesh, grid_t grid,  float *c, float mask, int algorithm, int verbose);

extern int fe_minmax(const mesh_t & mesh,frame_t & frame);

extern int fe_edgetable_T(mesh_t *mesh, int task, bool debug) ;
extern int fe_edgetable_Q(mesh_t *mesh, int verbose) ;
extern int fe_edgetable(mesh_t *mesh, int task=0, int verbose=0, bool debug=true);
extern int update_edge_table(mesh_t *mesh);

extern int report_edge_table(int interior,int boundary,int weird,int presumed,int nelts);
extern int init_edge_table(mesh_t * mesh, int verbose=0);

extern int fe_order(mesh_t & mesh, int n, double x, double y,int rotation);
extern int fe_order2(mesh_t & mesh, int n, int from, int rotation);
extern int fe_SortEdges(mesh_t & mesh, int n, int edge, int rotation);

extern int fe_oppositeneighbour(mesh_t mesh, int n, int m0);
extern int fe_anb(const mesh_t & mesh, int a, int b);
extern int fe_isedge(mesh_t & mesh, int a, int b);
extern int fe_quadnodes(mesh_t & mesh, int a, int b, int q[4]);

extern int fe_is_neighbour(mesh_t *mesh, int n, int m);

extern int fe_cleanvertices(mesh_t *mesh, vector<int> & watch, bool debug);
extern int fe_cleanvertices(mesh_t *mesh, bool debug);

extern int fe_cleanisolated(mesh_t & mesh, int maxpercent);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  mesh vertex editing

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern int fe_setvertex(mesh_t *mesh, double lon, double lat, float z,int *ngh, int nngh, int target);
extern int fe_addvertex(mesh_t *mesh, double lon, double lat, float z=NAN,int *ngh=0, int nngh=0); 
extern int fe_AddVertexOnEdge(mesh_t *mesh, int edge);

extern int fe_addvertex3D(mesh_t *mesh, double lon, double lat, double h, double *zlevel,int *ngh, int nngh);

extern int fe_insertvertex(mesh_t *mesh, int n);
extern int fe_removevertex(mesh_t *mesh, int n);
extern int fe_displacevertex(mesh_t & mesh, int target, double lon, double lat);

extern int fe_mergevertex(mesh_t *mesh, int n, int m, int remove=1); //TL

extern int fe_disconnectvertices(mesh_t & mesh, int n1, int n2);
extern int fe_disconnectvertex  (mesh_t & mesh, int n1, int n2);
extern int fe_disconnectvertex  (mesh_t & mesh, int n);

extern int fe_connectvertices(mesh_t & mesh, int n1, int n2);

extern int fe_substitute_vertex(mesh_t & mesh, int n1, int n2, int target);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  mesh geometry

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern double fe_distance(const mesh_t & mesh, int n, int m);
extern double fe_distance(const mesh_t & mesh, int n, double, double);
extern double fe_distance(const discretisation_t & descriptor, int n, double t, double p);

extern double fe_angle(mesh_t & mesh, int m1, int m2, int m3);
extern double fe_angle_NoRecale(mesh_t & mesh, int m1, int m2, int m3);

extern double fe_angle_cartesian(mesh_t & mesh, int m1, int m2, int m3);
extern double fe_angle_cartesian(mesh_t & mesh, int m1, int m2, int m3, double *x, double *y);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  mesh topology

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern int fe_codetable_obsolete(mesh_t *mesh, int RecycleCodes, int StopOnPinched, int verbose);
extern int fe_codetable1(mesh_t *mesh, int RecycleCodes, int stopon_EdgeError, int stopon_PinchError);
extern int fe_codetable2(mesh_t *mesh, int RecycleCodes, int stopon_EdgeError, int stopon_PinchError, int verbose=0);
extern int fe_codetable(mesh_t & mesh, int RecycleCodes, int SetLimits, int verbose);

extern int fe_vertex_element_tables(mesh_t *mesh);
extern int fe_vertex_Etable(mesh_t & mesh, int n);

extern int fe_setOBCflags(mesh_t & mesh, int n1, int n2, int flag);

extern mesh_t fe_refine(mesh_t & mesh,int *selected, int limits_locked);

extern int smooth_density_clx(const grid_t & grid, float *density, float mask, float slopemax);
extern int smooth_density_new(const grid_t & grid, float *density, float mask, float slopemax);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  mesh merging

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern mesh_t fe_merge(const mesh_t work[2], double dmax, const char *filename=0, bool limited=false, bool debug=false);
extern mesh_t fe_merge(mesh_t & small, mesh_t & big, double dmax);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  mesh splitting

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern mesh_t *fe_split(mesh_t & mesh,int *selected,int *targeted, int *frontier, int size, string rootname="anonymous", bool debug=false);
extern mesh_t *fe_split(mesh_t & mesh, vector<plg_t> & polygons, bool exact, bool debug);

extern int fe_MeshCut(mesh_t & mesh, vector<plg_t> & polygons, bool exact, bool debug);
extern int fe_MeshCut(mesh_t & mesh, const char *polygonfile, bool exact, bool debug);

extern int fe_presplit(mesh_t & mesh, const char *polyggonfile, double threshold, string rootname);
extern int fe_presplit(mesh_t & mesh, vector<plg_t> & polygons, double threshold, string rootname);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  mesh reshapping

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern int fe_reshape(mesh_t & mesh, int npass, vector<int> vtargets);

extern int fe_reshapeall01(mesh_t mesh,int npass);
extern int fe_reshapeall02(mesh_t mesh,int npass);
extern int fe_reshapeall(mesh_t & mesh,int npass);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  mesh editing

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern int fe_cleave(mesh_t *mesh, int ndel, bool debug);
extern int fe_cleavelimit(mesh_t *mesh, int ndel);
extern int fe_cleave(mesh_t *mesh, int ndel,int & n1,int & n2, bool debug);

extern int fe_exchange(mesh_t & mesh, int e, bool reshape, bool debug);

extern int fe_connex(mesh_t &);

extern int fe_reducebw(mesh_t, int);
extern int fe_reducebw(mesh_t & mesh, int incr, int max_percent, const char *output);

extern int fe_GetConnex(mesh_t & mesh, vector <vector<int> > & partition);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  mesh construction

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

extern mesh_t fe_doubleT(const mesh_t & mesh);
extern mesh_t fe_doubleQ(const mesh_t & mesh);
extern mesh_t fe_double(const mesh_t & mesh, int);

extern mesh_t *fe_node2mesh(mesh_t mesh, bool debug);

extern mesh_t *fe_createnodes(plg_t *polygones, int npolygones, point_t *interiors, int ninteriors);
extern int fe_createnodes(const vector<plg_t> polygones, point_t *interiors, int ninteriors, mesh_t & mesh);

extern int fe_createnodes(vector<plg_t> & polygons, double *x, double *y, double *z, int ninteriors, mesh_t & mesh);
extern int fe_createnodes(plg_t *polygones, int npolygones, double *x, double *y, int ninteriors, mesh_t & mesh);

extern int fe_limit_interior (mesh_t,limit_t,double,double,int *);

extern int fe_triangulateio_init(triangulateio & out);
extern int fe_node2triangle(mesh_t &, triangulateio *);
extern int fe_triangle2mesh(triangulateio &, mesh_t *);
extern mesh_t *fe_triangulate(mesh_t & nodes, projPJ projection, int verbose, bool debug);
extern int fe_triangulate(mesh_t & nodes, projPJ projection, mesh_t & finished, int verbose, bool debug);
extern void fe_FreeTriangleIO(triangulateio & io);

extern grid_t get_grid_n(const mesh_t & mesh);
extern grid_t get_grid_d(const mesh_t & mesh);

extern int fe_test_area(const mesh_t & mesh,int m);

extern int fe_intpl_LGP1(mesh_t & mesh, float *buffer, float mask, double t, double p, int m, float *z);
extern int fe_intpl_LGP1(mesh_t & mesh, double *buffer, double mask, double t, double p, int m, double *z);
extern int fe_intpl_LGP1(mesh_t & mesh, fcomplex *buffer, fcomplex mask, double t, double p, int m, fcomplex *z);

extern int fe_intpl_QLP1(const mesh_t & mesh,const float *buffer,float mask,double t,double p,int m,float *z);
extern int fe_intpl_QLP1(const mesh_t & mesh,const double *buffer,double mask,double t,double p,int m,double *z);
extern int fe_intpl_QLP1(const mesh_t & mesh,const complex<float> *buffer,complex<float> mask,double t,double p,int m,complex<float> *z);
extern int fe_intpl_QLP1(const mesh_t & mesh,const complex<double> *buffer,complex<double> mask,double t,double p,int m,complex<double> *z);

extern int fe_intpl_CQN1(const mesh_t & mesh,const float  *buffer,          float mask,           double t, double p, int m, float  *z);
extern int fe_intpl_CQN1(const mesh_t & mesh,const double *buffer,          double mask,          double t, double p, int m, double *z);
extern int fe_intpl_CQN1(const mesh_t & mesh,const complex<float> *buffer,  complex<float> mask,  double t, double p, int m, complex<float> *z);
extern int fe_intpl_CQN1(const mesh_t & mesh,const complex<double> *buffer, complex<double> mask, double t, double p, int m, complex<double> *z);

// make default mask value NAN (inactive)
template<typename T> int fe_intpl_LGP1(mesh_t mesh, T *buffer, double t, double p, int m, T *z)
{return fe_intpl_LGP1(mesh,buffer,NAN,t,p,m,z);}

extern int fe_interpolate2D(const mesh_t & mesh,int discretisation,const float *buffer,float mask,double t,double p,const int m,float *z);
extern int fe_interpolate2D(const mesh_t & mesh,int discretisation,const double *buffer,double mask,double t,double p,const int m,double *z);
extern int fe_interpolate2D(const mesh_t & mesh,int discretisation,const complex<float> *buffer,complex<float> mask,double t,double p,const int m,complex<float> *z);
extern int fe_interpolate2D(const mesh_t & mesh,int discretisation,const double *buffer,double mask,double t,double p,const int m,complex<double> *z);
extern int fe_interpolate2D(const mesh_t & mesh,int discretisation,const complex<double> *buffer,complex<double> mask,double t,double p,const int m,complex<double> *z);

extern double fe_extrapolate(const discretisation_t &descriptor,const complex<float> *tide,double x,double y,double dmax,complex<float> *zz,int verbose);
extern double fe_extrapolate(const discretisation_t &descriptor,const complex<double> *tide,double x,double y,double dmax,complex<double> *zz,int verbose);

extern int fe_sliceH(const mesh_t & mesh, int discretisation, double *depths, double level, float *buffer, float mask, float *slice);
extern int fe_sliceH(const mesh_t & mesh, int discretisation, double *depths, double level, double *buffer, double mask, double *slice);
extern int fe_sliceH(const mesh_t & mesh, int discretisation, double *depths, double level, double *buffer, double mask, double *slice);
extern int fe_sliceH(const mesh_t & mesh, int discretisation, double *depths, double level, complex<double> *buffer,complex<double> mask,complex<double> *slice);

extern int fe_sliceV(const mesh_t&mesh,const grid_t&grid,const list_t&list,int discretisation,double*depths,double*x,double*y,int nloc,float*buffer,float mask,grid_t*slice_grid,float**slice,int verbose);
extern int fe_sliceV(const mesh_t&mesh,const grid_t&grid,const list_t&list,int discretisation,double*depths,double*x,double*y,int nloc,double*buffer,double mask,grid_t*slice_grid,double**slice,int verbose);
extern int fe_sliceV(const mesh_t&mesh,const grid_t&grid,const list_t&list,int discretisation,double*depths,double*x,double*y,int nloc,double*buffer,double mask,grid_t*slice_grid,complex<double>**slice,int verbose);
extern int fe_sliceV(const mesh_t&mesh,const grid_t&grid,const list_t&list,int discretisation,double*depths,double*x,double*y,int nloc,complex<double>*buffer,complex<double> mask,grid_t*slice_grid,complex<double>**slice,int verbose);


extern void fe_set_triangles(mesh_t *mesh,const int *corners,int index0=0);
extern void fe_set_triangles(mesh_t *mesh,const double *corners,int index0=0);

extern void fe_set_quadrangles(mesh_t *mesh,const int *corners,int index0=0);
extern void fe_set_quadrangles(mesh_t *mesh,const double *corners,int index0=0);

extern int fe_readgeometry(const char *filename,mesh_t *mesh,int verbose=1) ;
extern void fe_init_from_elements(mesh_t *mesh,int verbose=0);
extern int fe_readmesh3d(const char *filename,mesh_t *mesh, int verbose) ;

extern int fe_savemeshNC2D(const char *filename,mesh_t & mesh,int option);
extern int fe_savemeshNC3D(const char *filename,mesh_t & mesh,int option);

extern int fe_savediscretisation(const char *, mesh_t &, discretisation_t &, int verbose=1);

extern int fe_savemesh3d(const char *filename, mesh_t & mesh, int option);
extern int fe_savemesh3d(const char *filename, mesh_t & mesh, int option, bool edges, bool levels, bool code, bool bathymetry);


extern int fe_CheckDimensions(const char *filename, const mesh_t & mesh);


// extern int archiving_UGdummy2D(const char *localname, mesh_t mesh,const char *name, double *buffer, int frame, int discretisation, int verbose=1);
// extern int archiving_UGdummy2D(const char *localname, mesh_t mesh,const char *name, double *buffer, int discretisation);
// extern int archiving_UGdummy2D(const char *localname, mesh_t mesh, const char *name1,const char *name2, complex<double> *buffer, int frame, int discretisation);
// extern int archiving_UGdummy2D(const char *localname, mesh_t mesh, const char *name1,const char *name2, complex<float> *buffer, int frame, int discretisation);
// extern int archiving_UGdummy2D(const char *localname, mesh_t mesh, const char *name1,const char *name2, complex<double> *buffer, int discretisation);
// extern int archiving_UGdummy2D(const char *localname, mesh_t mesh, const char *name1,const char *name2, complex<float> *buffer, int discretisation);

extern int archiving_UGdummy2D(const char *, const mesh_t &, const char *, const char *, const float *,  float,  int);
extern int archiving_UGdummy2D(const char *, const mesh_t &, const char *, const char *, const double *, double, int);
extern int archiving_UGdummy2D(const char *, const mesh_t &, const char *, const char *, const float *,  float,  int, int);
extern int archiving_UGdummy2D(const char *, const mesh_t &, const char *, const char *, const double *, double, int, int);

extern int archiving_UGdummy2D(const char *, const mesh_t &, const char *, const char *, const char *, complex <float> *, complex <float>, int discretisation);
extern int archiving_UGdummy2D(const char *, const mesh_t &, const char *, const char *, const char *, complex<double> *, complex<double>, int);
extern int archiving_UGdummy2D(const char *, const mesh_t &, const char *, const char *, const char *, complex<double> *, complex<double>, int, int);

extern int read_UG2Datlas(char *filename,mesh_t & mesh,fcomplex * &tide,fcomplex & mask,char *discretisation,int iteration,int verbose);

extern int archiving_UGdummy3D(const char *, mesh_t &, const char *, double **, int, int, int, int);
extern int archiving_UGdummy3D(const char *, mesh_t &, const char *, const char *, double **, int, int);
extern int archiving_UGdummy3D(const char *, mesh_t &, const char *, const char *, complex <double> **, int, int, int);
extern int archiving_UGdummy3D(const char *, mesh_t &, const char *, double **, int, int,int);


extern mesh_t *fe_partition_01(mesh_t mesh,int *partition,int npartition);
extern mesh_t *fe_partition_02(mesh_t mesh,int *partition,int npartition);
extern mesh_t *fe_partition_03(mesh_t mesh,int *partition,int npartition);
extern int *fe_read_npartition(const char *filename, mesh_t mesh);


extern int build_vindex(mesh_t mesh);

// class connexion_t {
// private :
// public :
//   int **incidence, *cardinal;
//   connexion_t() {
//     incidence=0;
//     cardinal=0;
//     }
//
// };
extern int fe_readdiscretisation(const char *filename,mesh_t *mesh, int fmt, int discretisation);

// from discretisation.cpp
extern const discretisation_t *get_descriptor_address(const mesh_t &mesh,int discretisation);
extern discretisation_t *get_descriptor_address(mesh_t &mesh,int discretisation);
extern const discretisation_t & get_descriptor(const mesh_t & mesh, int discretisation);
extern int paire_discretisation_init(std::string discretisation);
extern int paire_discretisation_id(int paire, int *z_discretisation, int *u_discretisation);
extern int paire_definition(int z_discretisation,int u_discretisation);
extern int paire_dimension(const mesh_t& mesh, int paire, int *zdim, int *udim);
extern int vector_dimension(mesh_t mesh,int discretisation, int *dim);
extern int discretisation_from_name(const char *name,int verbose=1);

// from discretisation-initialise.cpp
extern int paire_discretisation(const mesh_t & mesh, int paire, discretisation_t *z_descriptor,discretisation_t *u_descriptor);

/* *-----------------------------------------------------------------------------
Parallel computing routines */

extern int fe_crosstables(mesh_t *);

extern int fe_edge_crosstables01(mesh_t *);
 
extern int fe_vertex_crosstables01(mesh_t *);
extern int fe_vertex_crosstables02(mesh_t *);
extern int fe_element_surrounding( mesh_t);

extern int fe_Vertex_EdgesXtables(mesh_t *mesh);


extern int fe_element_crosstables(mesh_t);

/* *-----------------------------------------------------------------------------
Computational discretisation routines */

extern discretisation_t init_discretisation(const mesh_t &mesh, int type);
extern  int discretisation_init(mesh_t* mesh, int discretisation, int massmatrix = 1, int context=0);

extern int discretisation_PaireIdFromName(const string & name);
extern const char *discretisation_name(int discretisation);

/* *-----------------------------------------------------------------------------
  Packed matrix index initialisation */
extern int init_packed_LGP0xLGP0_0(mesh_t, ordering_t *);
extern int init_packed_LGP0xLGP0_1(mesh_t, ordering_t *);
extern int init_packed_LGP0xLGP0_2(mesh_t, ordering_t *);

extern int init_packed_LGP1xLGP0(mesh_t , ordering_t *);
extern int init_packed_LGP1xNCP1(mesh_t, ordering_t *);
extern int init_packed_LGP1xLGP2(mesh_t , ordering_t *);

extern int init_packed_LGP1xLGP1_0(mesh_t, ordering_t *);
extern int init_packed_LGP1xLGP1_1(mesh_t, ordering_t *);
extern int init_packed_LGP1xLGP1_2(mesh_t, ordering_t *);
extern int init_packed_LGP1xLGP1_3(mesh_t, ordering_t *,paire_t *, int);

extern int init_packed_NCP1xLGP0(mesh_t, ordering_t *);
extern int init_packed_NCP1xLGP1(mesh_t, ordering_t *);
extern int init_packed_NCP1xLGP2(mesh_t, ordering_t *);

extern int init_packed_NCP1xLGP0_02(mesh_t, ordering_t *, int);
extern int init_packed_NCP1xNCP1_02(mesh_t, ordering_t *, int);

extern int init_packed_LGP0xNCP1(mesh_t, ordering_t *);
extern int init_packed_LGP0xLGP1(mesh_t, ordering_t *);
extern int init_packed_LGP0xLGP0_LGP1(mesh_t mesh, ordering_t *ordering);

extern int init_packed_NCP1xNCP1(mesh_t mesh, ordering_t *ordering);

extern int init_packed_LGP2xLGP1 (mesh_t, ordering_t *);
extern int init_packed_LGP2xDGP1 (mesh_t, ordering_t *);
extern int init_packed_LGP2xNCP1 (mesh_t, ordering_t *);

extern int init_packed_LGP2xLGP2_0(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);
extern int init_packed_LGP2xLGP2_1(mesh_t mesh, ordering_t *ordering);
extern int init_packed_LGP2xLGP2_2(mesh_t mesh, ordering_t *ordering);

extern int init_packed_DGP1xLGP2(mesh_t mesh, ordering_t *ordering);

extern int init_packed_CQP0xCQP0(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);
extern int init_packed_CQP1xCQP1(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);
extern int init_packed_CQP1xCQP0(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);
extern int init_packed_CQN1xCQP0(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);
extern int init_packed_CQP0xCQP0_01(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);
extern int init_packed_CQP0xCQN1(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);
extern int init_packed_CQN1xCQN1(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);

extern int initialise_MassMatrix(mesh_t & mesh, discretisation_t & descriptor, hypermatrix_t *matrix, int verbose);
extern int dMassMatrix(mesh_t & mesh, discretisation_t & descriptor, hypermatrix_t *matrix);

/* *-----------------------------------------------------------------------------
  mesh generation and verification */


extern int *fe_read_epartition(const char *filename, const mesh_t & mesh);

extern int fe_chkvertex_consistency(mesh_t *mesh, int force);
extern int fe_chkvertex_00(mesh_t mesh, int *selected, int *histogram);
extern int fe_chkvertex_01(mesh_t mesh, int *selected, int *histogram);

extern int fe_Echk_crossing (mesh_t mesh,int *targeted);
extern int fe_Echk_angles(mesh_t & mesh, int *selected, double threshold, int mode, bool exterior_only, bool debug);

extern int fe_Tchk_angles(mesh_t & mesh, int *selected, double threshold, bool debug);

extern int fe_Tchk_AspectRatio(mesh_t & mesh, int *selected, double min_ratio,float *flags);
extern double fe_AspectRatio(mesh_t mesh, int m);
extern int fe_chkelement_02(mesh_t mesh, int *selected);

extern int fe_selectedges_01(const mesh_t & mesh, const char *polygonfile, int *selected, bool polar=false);
extern int fe_selectedges_01(const mesh_t & mesh, vector<plg_t> &polygons, int *selected, bool polar=false);
extern int fe_selectedges_01(const mesh_t & mesh, vector<plg_t> & polygons, int *selected, int option, bool polar);

extern int fe_selectedges_02(const mesh_t & mesh, int *selected, bool);
extern int fe_selectedges_03(const mesh_t & mesh, int *selected, bool);
extern int fe_selectedges_04(const mesh_t & mesh, double dmax, int *selected);

extern int fe_CureOverconnections(mesh_t *mesh, int, int, vector<int> & watch, bool);
extern int fe_fix_RiverMouth(mesh_t *mesh);
extern int fe_clean03(mesh_t *mesh, double);

extern  int fe_supress_diamond(mesh_t *mesh, int targetted, int nloopmax, bool debug);


extern int fe_load_criteria(const char *filename,criteria_t *criteria, int verbose=0);
extern int fe_save_criteria(const char *filename,criteria_t criteria);

extern int fe_InitCriteria(const char *specification,criteria_t & criteria);

extern int fe_loadboundaries(char *filename, int format, plg_t **polygones, int *npolygones);

extern int fe_limits2poly(const mesh_t & mesh, plg_t **polygones, int *npolygones, bool debug);
extern int fe_limits2poly(const mesh_t & mesh, vector<plg_t> & polygons, char *flag=0, bool debug=false);

extern int defproj(grid_t *grid, plg_t *polygones, int npolygones);
extern mesh_t fe_nodit(criteria_t criteria, grid_t sgrid, float *density, float mask, plg_t *polygones, int npolygones, bool debug, int verbose=1);

extern int initialise_meshgrid(const char *output,const criteria_t & criteria, grid_t *grid, float **density, float *mask, plg_t *polygones, int npolygones);

extern int fe_ComputeMeshsize (const char *output, const criteria_t & criteria, grid_t *grid, float **density, float *mask, plg_t *polygones, int npolygones, int initialize, bool debug);
extern int fe_ReloadMeshsize  (const char *filename, grid_t *grid, float **density, float *mask, plg_t *polygones, int npolygones);
extern int fe_ModifyMeshsize  (const char *output, criteria_t criteria, grid_t grid, float *prior, float mask, plg_t *polygones, int npolygones);

extern int fe_DecimateInterior(criteria_t & criteria, mesh_t & mesh, grid_t & grid, float *density, float mask, plg_t *polygons, int npolygons);

extern int fe_defgrid(grid_t *grid, vector<plg_t> & polygons, vector<plg_t> & limits, double step);
extern int fe_defgrid(grid_t *grid, plg_t *polygons, int npolygons, vector<plg_t> & limits, double step);

extern int mesh_resolution(mesh_t & mesh, float* & resolution_LGP0);
extern int mesh_resolution(mesh_t & mesh);
extern float *mesh_LGP0resolution(mesh_t & mesh);

extern int fe_ZGrid2Quadrangle(const grid_t & zgrid, const grid_t & fgrid, char  *landmask, float *topo, float topomask, mesh_t & mesh, int FixPinched, int verbose);
extern int fe_FGrid2Quadrangle(grid_t grid, grid_t cgrid, char*,  float *, float, mesh_t & mesh, const char*echo=0, bool debug=true);
extern int quadrangle_list (mesh_t & mesh, int Minsize=0, int MaxSize=6);
extern int quadrangle_ImportStructured(const char *gridfile, const char *maskfile, mesh_t & mesh, metagrid_t & input, double tag, string rootname, const char *wprojection, bool debug);
extern int quadrangle_ImportStructured(const char * filename, mesh_t & mesh, metagrid_t meta, Cgrid_t & Cgrid, double tag);

extern int fe_mimicgrid(grid_t vgrid, mesh_t *mesh, float *v_landmask, float *v_topo, grid_t zgrid, float *z_landmask, const char *filename, bool debug);
extern int fe_mimicgrid(grid_t vgrid, mesh_t *mesh, char  *v_landmask, float *v_topo, grid_t zgrid, char  *z_landmask, const char *filename, bool debug);

// extern int projection_LGP1xLGP0(double *in,double *out,mesh_t & mesh);
// extern int projection_LGP1xLGP0(complex<double> *in,complex<double> *out,mesh_t & mesh);

extern int projection_NCP1xLGP1(double *,          double *,          mesh_t  &);
extern int projection_NCP1xLGP1(complex<double> *, complex<double> *, mesh_t  &);
extern int projection_NCP1xLGP0(double *, double *, mesh_t &);

extern int projection_LGP1xLGP0(double *, double *, mesh_t &);
extern int projection_LGP0xLGP1(double *, double *, mesh_t &);
extern int projection_LGP0xLGP1(complex<double> *in, complex<double> *out, mesh_t & mesh);

extern int projection_LGP1xLGP2(double *, double *, mesh_t &);
extern int projection_NCP1xLGP2(double *, double *, mesh_t &);
extern int projection_LGP0xLGP2(double *, double *, mesh_t &);

extern int projection_CQP1xCQP0(mesh_t & mesh, double *in, double *out);


extern int projection_LGP0(float *,           float *,           mesh_t &, int);
extern int projection_LGP0(double *,          double *,          mesh_t &, int);
extern int projection_LGP0(complex<double> *, complex<double> *, mesh_t &, int);

extern int projection_DGP1(double *, double *, mesh_t &, int);

extern int projection_DNP1(double *,          double *,          mesh_t &, int);
extern int projection_DNP1(complex<double> *, complex<double> *, mesh_t &, int);

extern int projection_LGP1(float *,           float *,          mesh_t &, int);
extern int projection_LGP1(double *,          double *,          mesh_t &, int);
extern int projection_LGP1(complex<double> *, complex<double> *, mesh_t &, int);

extern int projection_LGP2(float *,           float *,           mesh_t &, int);
extern int projection_LGP2(double *,          double *,          mesh_t &, int);
extern int projection_LGP2(complex<double> *, complex<double> *, mesh_t &, int);

extern int projection_NCP1(double *, double *, mesh_t &, int);

extern int fe_projection(mesh_t &, float *,           int, float *,           int);
extern int fe_projection(mesh_t &, double *,          int, double *,          int);
extern int fe_projection(mesh_t &, complex<float> *,  int, complex<float> *,  int);
extern int fe_projection(mesh_t &, complex<double> *, int, complex<double> *, int);

extern int fe_Mesh2DElasticityT(mesh_t & mesh, double* & K, bool debug);

extern int fe_initcodes_01(mesh_t & mesh, int & nboundaries, int verbose);
extern int fe_initcodes_02(mesh_t & mesh, int & nboundaries, int verbose);
extern int fe_initcodes(mesh_t & mesh, int & nboundaries, int verbose);

extern int fe_delaunay_create(const char *rootname, mesh_t & mesh, mesh_t & delaunay, bool meshout, bool debug);
extern int fe_delaunay_improve(const char *rootname, mesh_t & mesh, double ElementMaxSize, double ElementMaxRatio, bool debug);

extern int fe_cartesian(mesh_t & mesh, const char *projection);
extern int fe_spherical(mesh_t & mesh, const char *projection);


#endif
