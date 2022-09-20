#ifndef TUGO_PROTOTYPES_H
#define TUGO_PROTOTYPES_H

#include "config.h"

#include "tugo-classes.h"

#include "fe.h"

#include "cefmo.h"

/**---------------------------Function prototypes---------------------------*/

/**-----------------------------------------------------------------------------
  Computational discretisation */
extern int paire_dimension(const mesh_t & mesh,int paire, int *zdim,int *udim);

/**-----------------------------------------------------------------------------
  2D open boundary conditions */

extern int dBoundaryMatrix(mesh_t mesh, int, double *A, int *pivot, int hbw);
extern void build_boundary_conditions(mesh_t mesh, state2d_t * P1_state,parameter_t  P1_data);

extern int VelocityBCs_init(specification_obsolete_t spec);

extern int GetBdyObs(const char *, double, int verbose, int debug);
extern int TidalBoundaryConditions(mesh_t & mesh, discretisation_t &, specification_obsolete_t &);

extern void process_boundarycode(mesh_t &, int, bool);

extern int boundary_normal(specification_obsolete_t & spec, mesh_t & mesh, int node, double *nx, double *ny, double *L);
extern int boundary_normal(specification_obsolete_t & spec, mesh_t & mesh, int node);

extern void treat_boundarycode(mesh_t& mesh, specification_obsolete_t *spec, int discretisation, int code);


/**-----------------------------------------------------------------------------
  Linear algebra processing routines */

extern void VertexGradient_operator(const mesh_t & mesh);


/**-----------------------------------------------------------------------------
  Model configuration routines */
extern int configure_init(char *filename,int ,tugo_cfg_C *, int delta=0);

extern int hugo_terminate(void);

extern int allocate_action2D(mesh_t & mesh, int level);

/**-----------------------------------------------------------------------------
 Runtime diagnostics */

extern double SquareDst(double, double, double, double);

extern void update_uNodes_elevation(mesh_t &, state2d_t &, parameter_t &);


/**-----------------------------------------------------------------------------
  3D horizontal momentum diffusion */

extern void initialize_horizontaldiffusion(void);

/**-----------------------------------------------------------------------------
  Bottom friction */

extern void initialize_bottomfriction2D(mesh_t &, parameter_t &);

extern void initialize_u0 (mesh_t &mesh, parameter_t data2D, int, int);

/**-----------------------------------------------------------------------------
  Internal wave drag */
extern matrix2x2_t *SpectralDrag_init(discretisation_t, tidal_wave, double *, double *, parameter_t, double wdHmin, double wdHmax, int);

extern void initialize_IWdrag(mesh_t &, parameter_t &, int);


/**-----------------------------------------------------------------------------
  FE  non-standard processing routines */
extern int read_mesh(const char *type,const char *rootname,mesh_t & mesh);
extern int geometry(mesh_t *);

/**-----------------------------------------------------------------------------
  Non-standard IO routines */
extern void load_optionals(mesh_t & mesh, state2d_t & state2D, parameter_t & data2D,int fmt, int verbose);

extern void ShowDimension(FILE *);

/**-----------------------------------------------------------------------------
  Bottom topography */
extern int topo_init(mesh_t & mesh,parameter_t & IO_data);
extern void topo_gradient(mesh_t & mesh, parameter_t VertexData);

/**-----------------------------------------------------------------------------
  Tides and harmonic analysis */
extern void tidal_equilibrium_harmonics(tidal_wave w, parameter_t data,int nnodes, complex<double> *z);

extern int initialize_LSA (mesh_t &, int, parameter_t *, spectrum_t &, int);

extern int load_wave_list(char *);

/**-----------------------------------------------------------------------------
  Energy and mass diagnostics */

extern zmatrix2x2_t *spectral_friction01(mesh_t &, complex <double> *, complex <double> *, int);
extern zmatrix2x2_t *spectral_friction02(mesh_t &, complex <double> *, complex <double> *, int);

#   ifndef ADDWAVE
#      define  ADDWAVE
extern int addwave(spectrum_t *, tidal_wave, int);
#   endif
       // ADDWAVE



extern int projection_CQP1xCQP0(mesh_t & mesh, double *in, double *out);


extern  int projection_LGP1(double *,          double *,          mesh_t &, int);
extern  int projection_LGP1(complex<double> *, complex<double> *, mesh_t &, int);

extern int initialise_UGO2D(mesh_t & mesh);
extern int allocate_UGO2D(mesh_t &);

/* vvvvvvvvv STS vvvvvvvvvvv */

#include "tugo-commons.h"

extern complex<float>  cconvert(hconstant_t x, int k);
extern complex<double> zconvert(hconstant_t x, int k);  

extern int coriolis_setup(discretisation_t & descriptor);
extern int coriolis_setup(parameter_t & data, int nnodes);

#include "zapper.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static int paire_dimension(mesh_t & mesh,int *zdim, int *udim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=paire_dimension(mesh,gSequentialPaire2D,zdim,udim);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static int paire_dimension(mesh_t & mesh,int *zdim,int *udim, int *gdim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int g_discretisation, status;

   status=paire_dimension(mesh,gSequentialPaire2D,zdim,udim);
   if(g2D_discretisation==AUTO) {
     g_discretisation=gradient_natural_discretisation(u2D_discretisation);
     }
   if(g2D_discretisation==INTRINSIC) {
     g_discretisation=gradient_natural_discretisation(u2D_discretisation);
     }
   else {
     g_discretisation=g2D_discretisation;
     }
   status=vector_dimension(mesh,g_discretisation,gdim);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static int paire_order(int *zorder, int *uorder)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Return dynamics equation order

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{

/**----------------------------------------------------------------------
  set matrix dynamical coefficients*/
  switch (we_type) {
    case GWE_CLX:
/**----------------------------------------------------------------------
      Generalized wave equation, Lynch & Gray [1979] formulation*/
      *zorder=1;
      *uorder=1;
      break;

    case GWE_XXX:
/**----------------------------------------------------------------------
      Generalized wave equation, revisited Lynch & Gray [1979] formulation*/
      *zorder=2;
      *uorder=1;
      break;

/**----------------------------------------------------------------------
    Explicit wave equation*/
    case WE_EXPLICIT:
      *zorder=1;
      *uorder=1;
      break;

/**----------------------------------------------------------------------
    Implicit wave equation*/
    case WE_IMPLICIT:
      *zorder=1;
      *uorder=1;
      break;

    default:
      check_error(-1, "illegal mode", __LINE__, __FILE__, 1);
      break;
  }

  return (0);
}


#endif /* TUGO_PROTOTYPES_H */
