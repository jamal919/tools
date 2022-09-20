
/*******************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

*******************************************************************************/

#include <stdio.h>
#include <string.h>
#include <time.h>

#ifdef HAVE_MPI
#include <mpi.h>
#endif

#include "fe.h"
#include "matrix.h"

static int sproduct2D=0;

static int const gLinearSolverID=0,gSubLinearSolverID=0; /// HERE !!!

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_Factorize(hypermatrix_t *M, int level)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int solver_id=-1,status=-1;
  bool force_duplicate=false;
  int verbose=0;
  bool debug=false;
  
  if(M->distributor->context==-1) M->distributor->context=SEQUENTIAL_COMPUTING;

//   switch (level) {
//     case 0:
//       solver_id=gLinearSolverID;
//       break;
//     default:
//       solver_id=gSubLinearSolverID;
//       break;
//     }
    
  if(solver_id==-1) {
    status=ishere("PASTIX");
    if(status==1) solver_id=SOLVER_ID_PASTIX;
    }
    
  if(solver_id==-1) {
    status=ishere("UMFPACK");
    if(status==1) solver_id=SOLVER_ID_UMFPACK;
    }
    
  if(solver_id==-1) solver_id=SOLVER_ID_SpDOMESTIC;
  
  status = LinearSystem_initialise(M, (int *) 0, 0, solver_id, force_duplicate, verbose, debug);
  
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_LGP0(mesh_t mesh, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Note:

  31/10/2010

  This approach is somehow a waste of memory, as the mass matrix is
  diagonal for LGP0.

  It is used to standardize the multi-discretisation capabilities in
  T-UGOm.

  Mignt be revised later on.

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
/**----------------------------------------------------------------------

  Mass matrix for LGP0 integrale scalar product

----------------------------------------------------------------------*/
  int count, m, n;
  double *M;
  double C, area;

  discretisation_t descriptor=mesh.LGP0descriptor;

  count = ordering.pointer[descriptor.nnodes];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.ntriangles; m++) {
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    M[ordering.pointer[m]] += area * C;
    if(M[ordering.pointer[m]]==0.0)
      TRAP_ERR_EXIT(-1,"anomaly in mass matrix\n");
    }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_LGP0(mesh_t mesh, discretisation_t & descriptor, hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on DGP1 nodes*/
  if(M->ordering->cardinal == NULL)
    status = init_packed_LGP0xLGP0_0(mesh, M->ordering);

  M->packed = MassPMatrix_LGP0(mesh, *(M->ordering));
  M->discretisation = LGP0;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/

  M->distributor=&(descriptor.distributor);
  M->distributor->context=SEQUENTIAL_COMPUTING;
  
  sprintf(msg,"%s level=%d context=%d","LGP0 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);
  
  printf("LGP0 mass matrix factorisation status: %d\n", status);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_LGP1_NQ(mesh_t & mesh, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------

  Mass matrix for P1 nodal quadrature (NQP1)

----------------------------------------------------------------------*/
{
  int count, i, j, k, m, n;
  int n1, n2;
  int row;
  double *M, p[3], q[3];
  double area;

  count = ordering.pointer[mesh.nvtxs];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes

  Check : MANDATORY!!!

  Notes

  Bug found 18/07/2008

  should be corrected to use nodal quadrature instead of integrale

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */

  for(m = 0; m < mesh.ntriangles; m++) {
    area = mesh.triangles[m].TrueArea;
    if(isnan(area)==1) TRAP_ERR_EXIT(-1, "triangle geometric metrix not initialised\n");
    for(i = 0; i < 3; i++) {
      n1 = mesh.triangles[m].vertex[i];
      for(k = 0; k < 3; k++)
        q[k] = 0;
      q[i] = 1.;
      for(j = 0; j < 3; j++) {
        n2 = mesh.triangles[m].vertex[j];
        for(k = 0; k < 3; k++)
          p[k] = 0;
        p[j] = 1.;
        for(k = 0; k < ordering.cardinal[n1]; k++) {
          if(ordering.incidence[ordering.pointer[n1] + k] == n2) {
            row = k;
            break;
          }
        }
        M[ordering.pointer[n1] + row] += 2. * area * integraleP1xP1(p, q);
      }
    }
  }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_LGP1_EI(mesh_t mesh, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------

  Mass matrix for P1 integrale quadrature (EIP1)

----------------------------------------------------------------------*/
  int count, i, j, k, m, n;
  int n1, n2;
  int row;
  double *M, p[3],q[3],r[3];
  double area;

  count = ordering.pointer[mesh.nvtxs];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.ntriangles; m++) {
//    area = mesh.triangles[m].TrueArea;
    for(i = 0; i < 3; i++) {
      n = mesh.triangles[m].vertex[i];
      r[i]=mesh.vertices[n].c;
      }
    area = mesh.triangles[m].Area;
    for(i = 0; i < 3; i++) {
      n1 = mesh.triangles[m].vertex[i];
      for(k = 0; k < 3; k++)
        q[k] = 0;
      q[i] = 1.;
      for(j = 0; j < 3; j++) {
        n2 = mesh.triangles[m].vertex[j];
        for(k = 0; k < 3; k++)
          p[k] = 0;
        p[j] = 1.;
        for(k = 0; k < ordering.cardinal[n1]; k++) {
          if(ordering.incidence[ordering.pointer[n1] + k] == n2) {
            row = k;
            break;
            }
          }
        if(k==ordering.cardinal[n1]) {
          TRAP_ERR_EXIT(-1, "ill-posed matrix\n");
          }
//        M[pointer[n1] + row] += 2. * area * integraleP1xP1(p, q);
        M[ordering.pointer[n1] + row] += 2. * area * fe_integraleLGP1xLGP1xLGP1_2D(p,q,r);
        }
      }
    }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_LGP1(mesh_t & mesh, discretisation_t& descriptor, hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on LGP1 nodes*/
  if(M->ordering->cardinal == NULL) {
    status = init_packed_LGP1xLGP1_0(mesh, M->ordering);
    }
 /*----------------------------------------------------------------------
  create packed mass matrix*/
  switch (sproduct2D) {
    case (NQUAD):
/*------------------------------------------------------------------------------
      nodal quadrature matrix */
      M->packed  = MassPMatrix_LGP1_NQ(mesh, *(M->ordering));
      break;

    case (INTGL):
/*------------------------------------------------------------------------------
      true integrale matrix */
      M->packed  = MassPMatrix_LGP1_EI(mesh, *(M->ordering));
      break;

    default:
      TRAP_ERR_EXIT(-1, "illegal scalar product\n");
    }

  M->discretisation = LGP1;

 /*----------------------------------------------------------------------
  factorize  packed mass matrix*/
  M->distributor=&(descriptor.distributor);

  sprintf(msg,"%s level=%d context=%d","LGP1 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);

  if(status!=0) printf("LGP1 mass matrix factorisation status: %d\n", status);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_DGP1(mesh_t mesh, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------

  Mass matrix for discontinuous LGP1 integrale scalar product

----------------------------------------------------------------------*/
  int count, i, j, k, m, n;
  int n1, n2;
  int row;
  double *M, p[6],q[6];
  double C, area;

  discretisation_t descriptor=mesh.DGP1descriptor;

  count = ordering.pointer[descriptor.nnodes];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.ntriangles; m++) {
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n1 = descriptor.NIbE[m][i];
      for(k = 0; k < descriptor.nnpe; k++)
        q[k] = 0;
      q[i] = 1.;
      for(j = 0; j < descriptor.nnpe; j++) {
        n2 = descriptor.NIbE[m][j];
        for(k = 0; k < descriptor.nnpe; k++)
          p[k] = 0;
        p[j] = 1.;
        for(k = 0; k < ordering.cardinal[n1]; k++) {
          if(ordering.incidence[ordering.pointer[n1] + k] == n2) {
            row = k;
            break;
            }
          }
        if(k==ordering.cardinal[n1]) {
          TRAP_ERR_EXIT(-1, "ill-posed matrix\n");
          }
        M[ordering.pointer[n1] + row] += 2. * area * C * fe_integraleLGP1xLGP1_2D(q,p);
        }
      }
    }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_DGP1(mesh_t mesh, discretisation_t& descriptor, hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

  extern int init_packed_DGP1xDGP1_0(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on DGP1 nodes*/
  if(M->ordering->cardinal == NULL)
    status = init_packed_DGP1xDGP1_0(mesh, descriptor, M->ordering);

  M->packed = MassPMatrix_DGP1(mesh, *(M->ordering));
  M->discretisation = DGP1;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/

  M->distributor=&(descriptor.distributor);

  sprintf(msg,"%s level=%d context=%d","DGP1 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;


/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);
  printf("DGP1 mass matrix factorisation status: %d\n", status);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_NCP1(mesh_t mesh, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------

  Mass matrix for discontinuous NCP1 integrale scalar product

----------------------------------------------------------------------*/
  int count, i, j, k, m, n;
  int n1, n2;
  int row;
  double *M, p[6],q[6];
  double C, area;

  discretisation_t descriptor=mesh.NCP1descriptor;

  count = ordering.pointer[descriptor.nnodes];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.ntriangles; m++) {
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n1 = descriptor.NIbE[m][i];
      for(k = 0; k < descriptor.nnpe; k++)
        q[k] = 0;
      q[i] = 1.;
      for(j = 0; j < descriptor.nnpe; j++) {
        n2 = descriptor.NIbE[m][j];
        for(k = 0; k < descriptor.nnpe; k++)
          p[k] = 0;
        p[j] = 1.;
        for(k = 0; k < ordering.cardinal[n1]; k++) {
          if(ordering.incidence[ordering.pointer[n1] + k] == n2) {
            row = k;
            break;
            }
          }
        if(k==ordering.cardinal[n1]) {
          TRAP_ERR_EXIT(-1, "ill-posed matrix\n");
          }
        M[ordering.pointer[n1] + row] += 2. * area * C * fe_integraleNCP1xNCP1_2D(q,p);
        }
      }
    }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_NCP1(mesh_t mesh, discretisation_t& descriptor, hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

  extern int init_packed_NCP1xNCP1(mesh_t mesh, ordering_t *ordering);

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on DGP1 nodes*/
  if(M->ordering->cardinal == NULL)
    status = init_packed_NCP1xNCP1(mesh, M->ordering);

  M->packed = MassPMatrix_NCP1(mesh, *(M->ordering));
  M->discretisation = NCP1;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/

  M->distributor=&(descriptor.distributor);

  sprintf(msg,"%s level=%d context=%d","NCP1 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);

  printf("NCP1 mass matrix factorisation status: %d\n", status);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_DNP1(mesh_t mesh, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------

  Mass matrix for discontinuous NCP1 integrale scalar product

----------------------------------------------------------------------*/
  int count, i, j, k, m, n;
  int n1, n2;
  int col;
  double *M, p[6],q[6];
  double C, area;

  discretisation_t descriptor=mesh.DNP1descriptor;

  count = ordering.pointer[descriptor.nnodes];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.ntriangles; m++) {
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n1 = descriptor.NIbE[m][i];
      for(k = 0; k < descriptor.nnpe; k++)
        q[k] = 0;
      q[i] = 1.;
      for(j = 0; j < descriptor.nnpe; j++) {
        n2 = descriptor.NIbE[m][j];
        for(k = 0; k < descriptor.nnpe; k++)
          p[k] = 0;
        p[j] = 1.;
        col=matrix_relative_pos(&ordering,n1,n2);
        M[ordering.pointer[n1] + col] += 2. * area * C * fe_integraleNCP1xNCP1_2D(q,p);
        }
      }
    }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_DNP1(mesh_t mesh, discretisation_t& descriptor, hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

  extern int init_packed_DNP1xDNP1_0(mesh_t &mesh, discretisation_t &descriptor, ordering_t *ordering);

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on DGP1 nodes*/
  if(M->ordering->cardinal == NULL)
    status = init_packed_DNP1xDNP1_0(mesh, descriptor, M->ordering);

  M->packed = MassPMatrix_DNP1(mesh, *(M->ordering));
  M->discretisation = DNP1;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  M->distributor=&(descriptor.distributor);

  sprintf(msg,"%s level=%d context=%d","DNP1 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;


/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);

  printf("DNP1 mass matrix factorisation status: %d\n", status);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_LGP2(mesh_t & mesh, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Mass matrix for continuous LGP2 integrale scalar product

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int count, i, j, m, n;
  int n1, n2;
  int row;
  double *M;
  double C, area;

  discretisation_t descriptor=mesh.LGP2descriptor;

  count = ordering.pointer[descriptor.nnodes];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.ntriangles; m++) {
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n1 = descriptor.NIbE[m][i];
      for(j = 0; j < descriptor.nnpe; j++) {
        n2 = descriptor.NIbE[m][j];
        row=matrix_relative_pos(&ordering,n1,n2);
        if(row==-1) {
          TRAP_ERR_EXIT(-1, "ill-posed matrix\n");
          }
        M[ordering.pointer[n1] + row] += 2. * area * C * fe_LGP2xLGP2_2D[i][j];
        }
      }
    }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_LGP2(mesh_t mesh, discretisation_t& descriptor, hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on LGP2 nodes*/
  if(M->ordering->cardinal == NULL)
    status = init_packed_LGP2xLGP2_0(mesh, descriptor, M->ordering);

  M->packed = MassPMatrix_LGP2(mesh, *(M->ordering));
  M->discretisation = LGP2;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/

  M->distributor=&(descriptor.distributor);

  sprintf(msg,"%s level=%d context=%d","LGP2 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);

  printf("LGP2 mass matrix factorisation status: %d\n", status);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_DGP2(mesh_t mesh, discretisation_t descriptor, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------

  Mass matrix for discontinuous DGP2 integrale scalar product

----------------------------------------------------------------------*/
  int count, i, j, m, n;
  int n1, n2;
  int col;
  double *M;
  double C, area;

  count = ordering.pointer[descriptor.nnodes];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.ntriangles; m++) {
    C=mesh.triangles[m].cosinus;
    area = mesh.triangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n1 = descriptor.NIbE[m][i];
      for(j = 0; j < descriptor.nnpe; j++) {
        n2 = descriptor.NIbE[m][j];
        col=matrix_relative_pos(&ordering,n1,n2);
//        M[ordering.pointer[n1] + col] += 2. * area * C * fe_integraleLGP2xLGP2_2D(i,j);
        M[ordering.pointer[n1] + col] += 2. * area * C * fe_LGP2xLGP2_2D[i][j];
        }
      }
    }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_DGP2(mesh_t mesh, discretisation_t& descriptor, hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on LGP2 nodes*/
  if(M->ordering->cardinal == NULL)
    status = init_packed_LGP2xLGP2_0(mesh, descriptor, M->ordering);

  M->packed = MassPMatrix_LGP2(mesh, *(M->ordering));
  M->discretisation = DGP2;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/

  M->distributor=&(descriptor.distributor);

  sprintf(msg,"%s level=%d context=%d","DGP2 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);

  printf("DGP2 mass matrix factorisation status: %d\n", status);

  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_CQP0(mesh_t & mesh, discretisation_t & descriptor, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Mass matrix for LGP0 integrale scalar product

----------------------------------------------------------------------*/
  int count, m, n;
  double *M;
  double C, area;

  count = ordering.pointer[descriptor.nnodes];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.nquadrangles; m++) {
    C=mesh.quadrangles[m].cosinus;
    area = mesh.quadrangles[m].Area;
    M[ordering.pointer[m]] += area * C;
    }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_CQP0(mesh_t & mesh, discretisation_t & descriptor, hypermatrix_t *M, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on DGP1 nodes*/
  if(M->ordering->cardinal == NULL)
    status = init_packed_CQP0xCQP0(mesh, descriptor, M->ordering);

  M->packed = MassPMatrix_CQP0(mesh, descriptor, *(M->ordering));
  M->discretisation = CQP0;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/

  M->distributor=&(descriptor.distributor);
  M->distributor->context=SEQUENTIAL_COMPUTING;

  sprintf(msg,"%s level=%d context=%d","QCP0 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);

  if(status!=0) printf("%s : CQP0 mass matrix factorisation status=%d\n", __func__, status);
  
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_CQP1(mesh_t & mesh, discretisation_t & descriptor, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Mass matrix for QP1 integrale scalar product

----------------------------------------------------------------------*/
  int count, i, j, m, n;
  int n1, n2;
  int col;
  double *M;
  double C, area;

  count = ordering.pointer[descriptor.nnodes];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.nquadrangles; m++) {
    C    = mesh.quadrangles[m].cosinus;
    area = mesh.quadrangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n1 = descriptor.NIbE[m][i];
      for(j = 0; j < descriptor.nnpe; j++) {
        n2 = descriptor.NIbE[m][j];
        col=vpos(n2,&(ordering.incidence[ordering.pointer[n1]]),ordering.cardinal[n1]);
        if(col==-1) {
          TRAP_ERR_EXIT(-1, "ill-posed matrix\n");
          }
        M[ordering.pointer[n1] + col] += C * fe_integrale_QP1xQP1xQP1_2D(i,j,mesh.quadrangles[m].J);
        }
      }
    }
    
//   for(m = 0; m < mesh.nquadrangles; m++) {
//     C    = mesh.quadrangles[m].cosinus;
//     area = mesh.quadrangles[m].Area;
//     for(i = 0; i < descriptor.nnpe; i++) {
//       n1 = descriptor.NIbE[m][i];
//       col=pos(n1,&(ordering.incidence[ordering.pointer[n1]]),ordering.cardinal[n1]);
//       if(col==-1) {
//         TRAP_ERR_EXIT(-1, "ill-posed matrix\n");
//         }
// //      M[ordering.pointer[n1] + col] += fe_QP1xQP1_2D[i][i];
//       M[ordering.pointer[n1] + col] += 1;
//       }
//     }
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_CQP1(mesh_t & mesh, discretisation_t & descriptor, hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on DGP1 nodes*/
  if(M->ordering->cardinal == NULL)
    status = init_packed_CQP1xCQP1(mesh, descriptor, M->ordering);

  M->packed = MassPMatrix_CQP1(mesh, descriptor, *(M->ordering));
  M->discretisation = CQP1;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/

  M->distributor=&(descriptor.distributor);

  sprintf(msg,"%s level=%d context=%d","QCP1 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);

  printf("CQP1 mass matrix factorisation status: %d\n", status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double *MassPMatrix_CQN1(mesh_t & mesh, discretisation_t & descriptor, ordering_t ordering)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/* *----------------------------------------------------------------------

  Mass matrix for QP1 integrale scalar product

----------------------------------------------------------------------*/
  int count, i, m, n;
  int n1;
  double *M;
  double C, area;

  count = ordering.pointer[descriptor.nnodes];

/*------------------------------------------------------------------------------
  allocate packed matrix*/
  M = new double[count];

  for(n = 0; n < count; n++)
    M[n] = 0;

  for(m = 0; m < mesh.nquadrangles; m++) {
    C    = mesh.quadrangles[m].cosinus;
    area = mesh.quadrangles[m].Area;
    for(i = 0; i < descriptor.nnpe; i++) {
      n1 = descriptor.NIbE[m][i];
      M[ordering.pointer[n1]] += C * area /4.;
      }
    }
    
  return (M);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix_CQN1(mesh_t & mesh, discretisation_t & descriptor, hypermatrix_t *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
  char msg[256];

/*------------------------------------------------------------------------------
  mass matrix needed to project any field on DGP1 nodes*/
  if(M->ordering->cardinal == NULL)
    status = init_packed_CQN1xCQN1(mesh, descriptor, M->ordering);

  M->packed = MassPMatrix_CQN1(mesh, descriptor, *(M->ordering));
  M->discretisation = CQN1;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/

  M->distributor=&(descriptor.distributor);

  sprintf(msg,"%s level=%d context=%d","QCN1 mass matrix",mesh.level,M->distributor->context);
  M->name=(string) msg;

/**----------------------------------------------------------------------
  factorize  packed mass matrix*/
  status=dMassMatrix_Factorize(M, mesh.level);

  printf("CQN1 mass matrix factorisation status: %d\n", status);
  
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int initialise_MassMatrix(mesh_t & mesh, discretisation_t & descriptor, hypermatrix_t *matrix, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=-1;
/**----------------------------------------------------------------------
  create packed mass matrix*/

  switch (descriptor.type) {
    case LGP0:
      printf("\ninitialise LGP0 mass matrix\n");
      status = dMassMatrix_LGP0(mesh, descriptor, matrix);
      break;

    case LGP1:
/*------------------------------------------------------------------------------
      nodal quadrature matrix */
      printf("\ninitialise LGP1 mass matrix\n");
      status = dMassMatrix_LGP1(mesh, descriptor, matrix);
      break;

    case DGP1:
/*------------------------------------------------------------------------------
      nodal quadrature matrix */
      printf("\ninitialise DGP1 mass matrix\n");
      status = dMassMatrix_DGP1(mesh, descriptor, matrix);
      break;

    case NCP1:
/*------------------------------------------------------------------------------
      nodal quadrature matrix */
      printf("\ninitialise NCP1 mass matrix\n");
      status = dMassMatrix_NCP1(mesh, descriptor, matrix);
      break;

    case DNP1:
/*------------------------------------------------------------------------------
      nodal quadrature matrix */
      printf("\ninitialise DNP1 mass matrix\n");
      status = dMassMatrix_DNP1(mesh, descriptor, matrix);
      break;

    case LGP2:
/*------------------------------------------------------------------------------
      true integrale matrix */
      printf("\ninitialise LGP2 mass matrix\n");
      status = dMassMatrix_LGP2(mesh, descriptor, matrix);
      break;

    case DGP2:
/*------------------------------------------------------------------------------
      true integrale matrix */
      printf("\ninitialise DGP2 mass matrix\n");
      status = dMassMatrix_DGP2(mesh, descriptor, matrix);
      break;

    case CQP0:
      if(verbose==1) printf("\ninitialise CQP0 mass matrix\n");
      status = dMassMatrix_CQP0(mesh, descriptor, matrix, verbose);
      break;

    case CQP1:
      printf("\ninitialise CQP1 mass matrix\n");
      status = dMassMatrix_CQP1(mesh, descriptor, matrix);
      break;

    case CQN1:
      printf("\ninitialise CQN1 mass matrix\n");
      status = dMassMatrix_CQN1(mesh, descriptor, matrix);
      break;

    default:
      TRAP_ERR_EXIT(-1, "mass matrix initialization, illegal discretisation\n");
    }

  if(verbose==1) {
    printf("matrix   size=%.3f Mo\n",  (double) matrix->ordering->size()*64/1.e+06);
    printf("ordering size=%.3f Mo\n\n",(double) matrix->ordering->size()*32/1.e+06);
    }
  return (status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int dMassMatrix(mesh_t & mesh, discretisation_t & descriptor, hypermatrix_t *matrix)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status, verbose=0;
  status=initialise_MassMatrix(mesh, descriptor, matrix, verbose);
  return (status);
 
}
