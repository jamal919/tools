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
#ifndef MASS_PROTO
#define MASS_PROTO

#include "tools-structures.h"

extern int init_packed_LGP1xLGP1_0(mesh_t mesh, ordering_t *ordering);
extern double *MassPMatrix_LGP1_01(mesh_t mesh, ordering_t ordering);
extern double *MassPMatrix_LGP1_02(mesh_t mesh, ordering_t ordering);
extern double *Identitymatrix(mesh_t mesh, ordering_t ordering);
extern int MassMatrix_factorisation(mesh_t *mesh, double *packed,ordering_t ordering,triplet *Triplet, solver_t **Solver);
extern int MassMatrix_solve(mesh_t mesh, triplet *Triplet, solver_t *Solver, double *B);
extern int dMassMatrix_LGP1(mesh_t mesh, discretisation_t descriptor, hypermatrix_t *M, int sproduct2D);
extern void fe_NCP1_to_LGP1(double *in, double *out, mesh_t mesh);
extern void fe_LGP0_to_LGP1(double *in, double *out, mesh_t mesh);
extern double *MassPMatrix_LGP2(mesh_t mesh, ordering_t ordering);
extern int dMassMatrix_LGP2(mesh_t mesh, discretisation_t descriptor, hypermatrix_t *M);
extern double *MassPMatrix_DGP1(mesh_t mesh, ordering_t ordering);
extern int dMassMatrix_DGP1(mesh_t mesh, discretisation_t descriptor, hypermatrix_t *M);
extern int dMassMatrix(mesh_t & mesh, discretisation_t & descriptor, hypermatrix_t *matrix);


#endif


