
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

***************************************************************************/

#ifndef __CEFMO
#define __CEFMO

#include "tugo-prototypes.h"

#include <map>

#define CEFMO_ARCHIVES_MINIMAL      0
#define CEFMO_ARCHIVES_STANDARD     1
#define CEFMO_ARCHIVES_DEBUGGING    2
#define CEFMO_ARCHIVES_EXPLORATION  3

#define TCLOSURE_CONSTANT    0
#define TCLOSURE_YAMADA      1
#define TCLOSURE_GALPERIN    2
#define TCLOSURE_GASPAR      3
#define TCLOSURE_TKE         4
#define TCLOSURE_KOMEGA      5
#define TCLOSURE_HOMOGENEOUS 6

#define PML_ATTENUATION

class tide2D_t {
private :
public :
/*-----------------------------------------------------------------------------------------------------
  physical parameters*/
  double *h_z,*h_u;                       /* elevation nodes bathymetry, velocity nodes bathymetry          */
  complex <double> *z,*u,*v;              /* elevation, currents                                            */
  complex <double> *LSA,*potential;       /* Loading and self-attraction potential, astronomical potential  */
  complex <double> *ASP;                  /* atmospheric surface pressure                                   */

  tide2D_t() {
    h_z=0;
    h_u=0;
    u=0;
    v=0;
    z=0;
    LSA=0;
    potential=0;
    ASP=0;
    }

  template <typename T> void erase(T *v) {
    if(*v!=0) delete [] (*v);
    *v=0;
    }

  void destroy() {
    erase(&h_z);
    erase(&h_u);
    erase(&u);
    erase(&v);
    erase(&z);
    erase(&LSA);
    erase(&potential);
    erase(&ASP);
    }

  void duplicate() {//unused
    }
};

class tide3D_t {
private :
public :
/**----------------------------------------------------------------------------
  physical parameters*/
  double *h_z,*h_u;
  complex <double> **z,**u,**v,**w;      /* level displacement, layer velocities             */
  complex <double> *u_mean,*v_mean;      /* mean velocities                                  */
  complex <double> *LSA,*potential;      /* loading/self-attraction, astronomical potential  */
  double **rho, **rho__;                 /* density at z-node and u-node                     */
  int n_Unodes,n_Znodes,nlevels,nlayers; /*                                                  */

tide3D_t() {
    h_z=0;
    h_u=0;
    u=0;
    v=0;
    w=0;
    u_mean=v_mean=0;
    z=0;
    LSA=0;
    potential=0;
    rho=0;
    rho__=0;
    n_Unodes=n_Znodes=nlayers=nlevels=-1;
    }
    
  void allocate(int n1, int n2, int n3, int n4) {
    if(n1<=0) return;
    this->n_Unodes=n1;
    this->n_Znodes=n2;
    u_mean=new complex <double> [n_Unodes];
    v_mean=new complex <double> [n_Unodes];
    if((n2>0) && (n4>0)) {
      z  =new complex <double> *[n_Znodes];
      rho=new double *[n_Znodes];
      for(size_t n=0;n<n_Znodes;n++) {
        z[n]  =new complex <double> [n4];
        rho[n]=new double [n4];
        }
      this->nlevels=n4;
      }
    if(n3>0) {
      u=new complex <double> *[n_Unodes];
      v=new complex <double> *[n_Unodes];
      for(size_t n=0;n<n_Unodes;n++) {
        u[n]=new complex <double> [n3];
        v[n]=new complex <double> [n3];
        }
      this->nlayers=n3;
      }
    }
    
  template <typename T> void erase(T *v) {
    if(*v!=0) delete [] (*v);
    *v=0;
    }
  template <typename T> void erase(T *v, int nvalues) {
    if(*v!=0) {
      for(size_t n;n<nvalues;n++) {
	if(v[n]!=0) delete [] v[n];
        }
      delete [] (*v);
      }
    *v=0;
    }

  void destroy() {
    erase(&h_z);
    erase(&h_u);
    erase(&u_mean);
    erase(&v_mean);
    erase(&z,n_Znodes);
    erase(&rho,n_Znodes);
    erase(&u,n_Unodes);
    erase(&v,n_Unodes);
    erase(&rho,n_Unodes);
//     erase(&w);
    erase(&LSA);
    erase(&potential);
    }

  void duplicate() {//unused
    }
};

class actionZ_t {
private :
public :
  complex <double>  **prx, **pry;      /* pressure gradient                         */
  complex <double>  *tpx,   *tpy;      /* potential                                 */
  complex <double>  *wdx,   *wdy;      /* wave drag                                 */
  complex <double>  *Fx,     *Fy;      /* miscellaneous                             */

  actionZ_t() {
    prx=pry=0;
    tpx=tpy=0;
    wdx=wdy=0;
    Fx=Fy=0;
    }

  template <typename T> void erase(T *v) {
    if(*v!=0) delete [] (*v);
    *v=0;
    }

  void destroy() {
    erase(&prx);
    erase(&pry);
    erase(&tpx);
    erase(&tpy);
    erase(&wdx);
    erase(&wdy);
    erase(&Fx);
    erase(&Fy);
    }

};

class cefmo_t {
private :
public :
/*-----------------------------------------------------------------------------
  physical parameters*/
  zmatrix2x2_t    *FrictionMatrix;
  matrix2x2_t     *DragMatrix;
  tidal_wave dominants[2];
  int   ndominants;
/**----------------------------------------------------------------------------
  PDE discrete operators matrices */
  hyperzmatrix_t   SpMatrix;
  hyperzmatrix_t   InvSpMomentum;
  complex <double> **InvSpMu,**InvSpMv;
  complex <double> **SpMu,**SpMv;
  hypermatrix_t  D1,D2;
  hyperzmatrix_t M1,M2;
  hypermatrix_t  G1,G2;
  hyperzmatrix_t DIVx,DIVy;
  hyperzmatrix_t GRDx,GRDy;
  double *u_LumpedMW;
  int (*SpGradient2D)       (const mesh_t &, cefmo_t &, parameter_t &);
  int (*SpDivergence2D)     (mesh_t &, cefmo_t &, parameter_t &);
  int (*SpDivergence2D_RHS) (mesh_t & mesh, complex <double> *, complex <double> *, complex <double> *, int nnodes, complex <double> *, const int);
  int (*SpAdvection2D)      (mesh_t & mesh, double *h, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *);
  int (*packed_ZBASExZBASE) (mesh_t &, discretisation_t &, ordering_t *);
  int (*packed_UBASExZBASE) (mesh_t &, discretisation_t &, ordering_t *);
  int (*packed_ZBASExUBASE) (mesh_t &, discretisation_t &, ordering_t *);
/**----------------------------------------------------------------------------
  OBCs */
  int *clamped, nclamped;
/**----------------------------------------------------------------------------
  discretisation */
  int paire;
  discretisation_t *u_descriptor;
  discretisation_t *z_descriptor;
/**----------------------------------------------------------------------------
  State and forcing vectors */
  tide2D_t  state;
  tide3D_t  state3D;
  actionZ_t action;
  complex <double>  *NLRhs;
  string rootname;

  cefmo_t() {
    FrictionMatrix=0;
    DragMatrix=0;
    clamped=0;
    nclamped=-1;
    paire=-1;
    NLRhs=0;
    u_descriptor=0;
    z_descriptor=0;
    u_LumpedMW=0;
    SpMu=SpMv=0;
    InvSpMu=InvSpMv=0;
    rootname="spectral";
    }

  template <typename T> void erase(T  &v) {
    if(v!=0) delete [] (v);
    v=0;
    }

  void destroy() {
    erase(FrictionMatrix);
    erase(DragMatrix);
    erase(clamped);
    state.destroy();
    action.destroy();
    nclamped=-1;
    paire=-1;
    erase(NLRhs);
    u_descriptor=0;
    z_descriptor=0;
    rootname.clear();
    }
};

/* STS: cmp_str already in functions.h */

typedef std::map<const char *, tide2D_t, cmp_str> atlas2D_t;
typedef std::map<const char *, tide3D_t, cmp_str> atlas3D_t;

// class compound_t {
// private :
// public :
// /*-----------------------------------------------------------------------------------------------------
//   physical parameters*/
//   tidal_wave  wave;
//   tide2D_t    *state;
//   spectrum_t  generating;
//   vector<string> formula;
//   int id,n;
// 
//   compound_t() {
//     n=0;
// //    formula=string("");
//     }
// };

extern void initialize_OPENMP();
extern int cefmo_solver(mesh_t *mesh, parameter_t data);

extern void accelere(double rn0,double rn1,double *rn2);

extern  void momentum2d_DGP1xLGP2_SpInitialise(tidal_wave, cefmo_t &, parameter_t &,mesh_t  &);
extern  void momentum2d_DGP1xLGP2_SpInitialise_iterate(tidal_wave,cefmo_t &,mesh_t &);

extern  complex<double>  *WE2D_DGP1xLGP2_SpInitialise(mesh_t &,tidal_wave,cefmo_t &,parameter_t &, ordering_t *);
extern  void WE2D_DGP1xLGP2_SpInitialise_iterate(mesh_t &,tidal_wave,cefmo_t &,hyperzmatrix_t);

extern  int  WE2D_CQN1xCQP0_SpDivergence(mesh_t &, complex <double> *, complex <double> *, complex <double> *, int, complex <double> *, const int);

extern  int  WE2D_DGP1xLGP2_SpDivergence_SEQ(mesh_t mesh, complex <double> *z, complex <double> *u, complex <double> *v, int nnodes, complex <double> *divergence);
extern  int  WE2D_DGP1xLGP2_SpDivergence_OMP(mesh_t mesh, complex <double> *z, complex <double> *u, complex <double> *v, int nnodes, complex <double> *divergence);
extern  int  WE2D_DGP1xLGP2_SpDivergence(mesh_t mesh, complex <double> *z, complex <double> *u, complex <double> *v, int nnodes, complex <double> *divergence);

extern  int  WE2D_DNP1xLGP2_SpDivergence_SEQ(mesh_t & mesh, complex <double> *z, complex <double> *u, complex <double> *v, int nnodes, complex <double> *divergence, const int);
extern  int  WE2D_DNP1xLGP2_SpDivergence_OMP(mesh_t & mesh, complex <double> *z, complex <double> *u, complex <double> *v, int nnodes, complex <double> *divergence, const int);
extern  int  WE2D_DNP1xLGP2_SpDivergence(mesh_t & mesh, complex <double> *z, complex <double> *u, complex <double> *v, int nnodes, complex <double> *divergence, const int);

extern  void WE2D_DGP1xLGP2_SpRHS(mesh_t, cefmo_t, parameter_t, complex<double>  *);
extern  void momentum_DGP1xLGP2_SpSolver(mesh_t, cefmo_t, parameter_t);

extern void cefmo_frictionRHS_OMP(spectrum_t, spectrum_t, double, atlas2D_t, parameter_t, int, complex<double> **, complex<double> **);
extern void cefmo_frictionRHS_SEQ(spectrum_t, spectrum_t, double, atlas2D_t, parameter_t, int, complex<double> **, complex<double> **);
extern void cefmo_frictionRHS(spectrum_t, spectrum_t, double, atlas2D_t, parameter_t, int, complex<double> **, complex<double> **);

extern  int momentum2d_CQN1xCQP0_SpAdvection(mesh_t &, double *, complex <double> *, complex <double> *,complex <double> *, complex <double> *, complex <double> *, complex <double> *);

extern int  momentum2d_DGP1xLGP2_SpAdvection_SEQ(mesh_t, double *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *);
extern int  momentum2d_DGP1xLGP2_SpAdvection_OMP(mesh_t, double *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *);
extern int  momentum2d_DGP1xLGP2_SpAdvection(mesh_t, double *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *);

extern int  momentum2d_DNP1xLGP2_SpAdvection_SEQ(mesh_t &, double *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *);
extern int  momentum2d_DNP1xLGP2_SpAdvection_OMP(mesh_t &, double *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *);
extern int  momentum2d_DNP1xLGP2_SpAdvection(mesh_t &, double *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *, complex <double> *);

extern void spectral_friction02(mesh_t, zmatrix2x2_t *, complex <double> *, complex <double> *, int);
extern void friction_coefficient(zmatrix2x2_t FrictionMatrix, double Cd, double r, double h, complex<double> *r1, complex<double> *r2, complex<double> *r3, complex<double> *r4);
extern void friction_coefficient(zmatrix2x2_t FrictionMatrix, double Cd, double r, double u0, double h, double ratio, complex<double> *r1, complex<double> *r2, complex<double> *r3, complex<double> *r4);

extern void cefmo_diffusion3D(spectrum_t, double, atlas3D_t *, parameter_t &, mesh_t &, cefmo_t &, int);
extern void cefmo_diffusion3D_OpenMP(spectrum_t, double, atlas3D_t *, parameter_t &, mesh_t &, cefmo_t &, int);
extern int  cefmo_diffusion3D_OpenMP_new(spectrum_t prediction, double duration,atlas3D_t *atlas, parameter_t & data, mesh_t & mesh, cefmo_t & cefmo, int iteration);

extern int  cefmo_topography(mesh_t &, tide2D_t &, parameter_t &, int);

extern int cefmo2D_DGP1xLGP2(mesh_t, parameter_t, cefmo_t, atlas2D_t *atlas);
extern int cefmo2D_DNP1xLGP2(mesh_t, parameter_t, cefmo_t, atlas2D_t *atlas);
extern int cefmo2D_NCP1xLGP2(mesh_t, parameter_t, cefmo_t, atlas2D_t *atlas);
extern int cefmo2D_LGP0xLGP1(mesh_t, parameter_t, cefmo_t, atlas2D_t *atlas,int already_initialised=0);
extern int cefmo2D_LGP2xLGP2(mesh_t, parameter_t, cefmo_t, atlas2D_t *atlas);

extern int SpGradient2D_LGP0xLGP1  (const mesh_t & mesh, cefmo_t & cefmo, parameter_t & data);
extern int SpDivergence2D_LGP0xLGP1(mesh_t & mesh, cefmo_t & cefmo, parameter_t & data);

extern int SpGradient2D_CQN1xCQP0   (const mesh_t & mesh, cefmo_t & cefmo, parameter_t & data);
extern int SpDivergence2D_CQN1xCQP0 (mesh_t & mesh, cefmo_t & cefmo, parameter_t & data);

extern int cefmo3D_DGP1xLGP2(mesh_t, parameter_t, cefmo_t, atlas3D_t *, atlas2D_t);
extern int cefmo3D_LGP0xLGP1(mesh_t & mesh, parameter_t & data, cefmo_t & cefmo, atlas3D_t *atlas3D, atlas2D_t & atlas2D);

extern int decode_compound(string s, spectrum_t *generating, int **signs);
extern int decode_compound(string s, tidal_wave **generating, int **signs, int *ngenerating);


extern  int SpectralFriction_init(discretisation_t descriptor, tidal_wave wave,parameter_t data, int nndes);
extern  int friction_ratio(mesh_t & mesh, discretisation_t & u_descriptor, float *Cd, double *ratio);

/**----------------------------------------------------------------------------
 generic calls */
extern void SpMomentum2D_solve(mesh_t &, cefmo_t &, complex<double> *,complex<double> *);

extern void SpArchive_Put2D(const mesh_t & mesh,const cefmo_t & cefmo,const parameter_t & data,const tidal_wave & wave,int iteration);
extern int  SpArchive_Get2D(const char *filename, int, cefmo_t cefmo, tide2D_t state);

extern int  SpArchive_Get3D(const char *filename, int, cefmo_t cefmo, tide3D_t state);

extern int  cefmo_allocate(mesh_t & mesh, cefmo_t & cefmo, int already_initialised);
extern int  SpWE2D_Terminate(mesh_t mesh,cefmo_t cefmo);
extern int  cefmo2D(mesh_t mesh, parameter_t data, cefmo_t cefmo, atlas2D_t *atlas,int already_initialised);

extern void cefmo2D_implicitFBCs01(complex<double> **SpMu, complex<double> **SpMv);
extern void cefmo2D_implicitFBCs02(cefmo_t & cefmo);

extern int cefmo3D_compute_LGP0xLGP1(mesh_t &, parameter_t &, cefmo_t &, tidal_wave, int, double *, complex<double> *, complex<double> *, int );
extern void momentum3D_SpSolver(mesh_t& mesh, cefmo_t& model, parameter_t data);
extern void pressure3Dbaroclinic_LGP0xLGP1(mesh_t& mesh, cefmo_t& model, parameter_t data);
extern void momentum3D_LGP0xLGP1_SpInitialise(mesh_t & mesh, cefmo_t & model, parameter_t & data,tidal_wave wave);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> T **NZ_2_ZN(mesh_t & mesh, T **buffer, int h_discretisation, int v_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n,hdim,vdim;
  double **tmp;
  
  discretisation_t descriptor=get_descriptor(mesh,h_discretisation);
  hdim=descriptor.nnodes;
  
  switch (v_discretisation) {
    case LAYERS:
      vdim=mesh.nlayers;
      break;
    case LEVELS:
      vdim=mesh.nlevels;
      break;
    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }
/**----------------------------------------------------------------------------
  re-organize buffer order to save layer slices */    
  tmp=new double*[vdim];
  for(k=0;k<vdim;k++) {
    tmp[k]=new double[hdim];
    for (n=0;n<hdim;n++) {
      tmp[k][n]=buffer[n][k];
      }
    }
  return(tmp);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int ZN_2_NZ(mesh_t & mesh, T **buffer, T **out, int h_discretisation, int v_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n,hdim,vdim;
  double **tmp;
  
  discretisation_t descriptor=get_descriptor(mesh,h_discretisation);
  hdim=descriptor.nnodes;
  
  switch (v_discretisation) {
    case LAYERS:
      vdim=mesh.nlayers;
      break;
    case LEVELS:
      vdim=mesh.nlevels;
      break;
    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }
/**----------------------------------------------------------------------------
  re-organize buffer order */    
  for(k=0;k<vdim;k++) {
    for (n=0;n<hdim;n++) {
      out[n][k]=buffer[k][n];
      }
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int erase_ZN(mesh_t & mesh, T ***buffer, int v_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n,vdim;
  double **tmp=*buffer;
    
  switch (v_discretisation) {
    case LAYERS:
      vdim=mesh.nlayers;
      break;
    case LEVELS:
      vdim=mesh.nlevels;
      break;
    default:
      check_error(-1, "discretisation not implemented in archive procedure", __LINE__, __FILE__, 1);
      break;
    }
  for(k=0;k<vdim;k++) {
    delete [] tmp[k];
    }
  delete [] tmp;
  return (status);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int erase_NZ(mesh_t & mesh, T ***buffer, int h_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int k,n,hdim;
  double **tmp=*buffer;
    
  discretisation_t descriptor=get_descriptor(mesh,h_discretisation);
  hdim=descriptor.nnodes;
  
  for(k=0;k<hdim;k++) {
    delete [] tmp[k];
    }
  delete [] tmp;
  return (status);

}


#endif
