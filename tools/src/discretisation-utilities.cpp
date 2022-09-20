
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

#include "config.h"

#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <cmath>
#include <map>

#include "tools-structures.h"
#include "fe.h"
#include "mass.h"

extern void fe_print(mesh_t mesh, int m);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int discretisation_from_name(const char *name,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///get discretisation from its name
/**
\param name discretisation name
\returns the discretisation code or #UNSET if unrecognised
*/
/*----------------------------------------------------------------------------*/
{

// This macro avoids typing errors :
#define caseReturn(x) if(!strcmp(name,#x))return x
  caseReturn(CQP0);
  caseReturn(CQP1);
  caseReturn(CQN1);
  caseReturn(LGP0);
  caseReturn(LGP1);
  caseReturn(DGP1);
  caseReturn(NCP1);
  caseReturn(DNP1);
  caseReturn(LGP2);
  caseReturn(DGP2);
#undef caseReturn
  TRAP_ERR_RETURN(UNSET,verbose,"unknown discretisation %s\n",name);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

const char *discretisation_name(int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///give char name of discretisation
/*----------------------------------------------------------------------------*/
{
  switch (discretisation) {
// This macro avoids typing errors :
#define caseReturn(x) case x:return #x
    caseReturn(CQP0);
    caseReturn(CQP1);
    caseReturn(CQN1);
    caseReturn(LGP0);
    caseReturn(LGP1);
    caseReturn(DGP1);
    caseReturn(NCP1);
    caseReturn(DNP1);
    caseReturn(LGP2);
    caseReturn(DGP2);
#undef caseReturn
    default:
      TRAP_ERR_RETURN(NULL,1,"unknown discretisation %d\n",discretisation);
  }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int discretisation_PaireIdFromName(const string & name)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int paire_id;
  
  if(name=="DGP1xLGP2") {
    paire_id=DGP1xLGP2;
    }
  else if(name=="NCP1xLGP2") {
    paire_id=NCP1xLGP2;
    }
  else if(name=="DNP1xLGP2") {
    paire_id=DNP1xLGP2;
    }
  else if(name=="LGP0xLGP1") {
    paire_id=LGP0xLGP1;
    }
  else if(name=="NCP1xLGP0") {
    paire_id=NCP1xLGP0;
    }
  else if(name=="NCP1xLGP1") {
    paire_id=NCP1xLGP1;
    }
  else if(name=="LGP2xLGP2") {
    paire_id=LGP2xLGP2;
    paire_id=IPGxLGP2;
    }
  else if(name=="CQP1xCQP0") {
    paire_id=CQP1xCQP0;
    }
  else if(name=="CQN1xCQP0") {
    paire_id=CQN1xCQP0;
    }
  else {
    check_error(-1, "illegal paire designation...", __LINE__, __FILE__, 1);
    }
  return(paire_id);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int gradient_natural_discretisation(int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int g_discretisation;

  g_discretisation=-1;
  switch (discretisation) {
    case LGP0:
      g_discretisation=NCP1;
      break;

    case DGP1:
      g_discretisation=DNP1;
      break;

    case LGP1:
      g_discretisation=LGP0;
      break;

    case NCP1:
      g_discretisation=LGP0;
      break;

    case DNP1:
      g_discretisation=DNP1;
      break;

    case LGP2:
      g_discretisation=DGP1;
      break;

    default:
      check_error(-1, "illicit discretisation", __LINE__, __FILE__, 1);
      break;
   }
  return(g_discretisation);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int paire_discretisation_id(int paire, int *z_discretisation,int *u_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  switch(paire) {
    case NCP1xLGP0:
/**----------------------------------------------------------------------
      NCP1xLGP0 structure initialisation*/
      *z_discretisation=LGP0;
      *u_discretisation=NCP1;
      break;

    case NCP1xLGP1:
/**----------------------------------------------------------------------
      NCP1xLGP1 structure initialisation*/
      *z_discretisation=LGP1;
      *u_discretisation=NCP1;
      break;

    case NCP1xLGP2:
/**----------------------------------------------------------------------
      NCP1xLGP1 structure initialisation*/
      *z_discretisation=LGP2;
      *u_discretisation=NCP1;
      break;

    case LGP0xLGP1:
/**----------------------------------------------------------------------
      LGP1xLGP0 structure initialisation*/
      *z_discretisation=LGP1;
      *u_discretisation=LGP0;
      break;

    case LGP1xLGP0:
/**----------------------------------------------------------------------
      LGP1xLGP0 structure initialisation*/
      *z_discretisation=LGP0;
      *u_discretisation=LGP1;
      break;

    case LGP1xLGP1:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *z_discretisation=LGP1;
      *u_discretisation=LGP1;
      break;

    case LGP1xLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP2 structure initialisation*/
      *z_discretisation=LGP2;
      *u_discretisation=LGP1;
      break;

    case DGP1xLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *z_discretisation=LGP2;
      *u_discretisation=DGP1;
      break;

    case DNP1xLGP1:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *z_discretisation=LGP1;
      *u_discretisation=DNP1;
      break;

    case DNP1xLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *z_discretisation=LGP2;
      *u_discretisation=DNP1;
      break;

    case LGP2xLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *z_discretisation=LGP2;
      *u_discretisation=LGP2;
      break;

    case IPGxLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *z_discretisation=LGP2;
      *u_discretisation=IPG;
      break;

    case CQP1xCQP0:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *z_discretisation=CQP0;
      *u_discretisation=CQP1;
      break;

    case CQN1xCQP0:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *z_discretisation=CQP0;
      *u_discretisation=CQN1;
      break;

    default:
      check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
      break;
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int paire_definition(int z_discretisation,int u_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int paire2D;

  if(z_discretisation==-1) {
    check_error(1, "z discretization not initialized", __LINE__, __FILE__, 1);
    }

  if(u_discretisation==-1) {
    check_error(1, "u discretization not initialized", __LINE__, __FILE__, 1);
    }

  if(z_discretisation==u_discretisation) {
    char *msg;
    asprintf(&msg,"WARNING : same discretisation on speed and elevation (%s) IS UNSTABLE ON 4-LINKED-NODES : MAKE SURE YOUR MESH DOES NOT HAVE ANY.",discretisation_name(z_discretisation));
    check_error(1, msg, __LINE__, __FILE__, 0);
    free(msg);
    }

  paire2D=-1;
  switch (u_discretisation) {
    case LGP0:
      switch (z_discretisation) {
        case LGP1:
          paire2D=LGP0xLGP1;
          break;
        default:
          check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
          break;
       }
      break;

    case DGP1:
      switch (z_discretisation) {
        case LGP2:
          paire2D=DGP1xLGP2;
          break;
        default:
          check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
          break;
       }
      break;

    case LGP1:
      switch (z_discretisation) {
        case LGP0:
          paire2D=LGP1xLGP0;
          break;
        case LGP1:
          paire2D=LGP1xLGP1;
          break;
        case LGP2:
          paire2D=LGP1xLGP2;
          break;
        default:
          check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
          break;
       }
      break;

    case NCP1:
      switch (z_discretisation) {
        case LGP0:
          paire2D=NCP1xLGP0;
          break;
        case LGP1:
          paire2D=NCP1xLGP1;
          break;
        case LGP2:
          paire2D=NCP1xLGP2;
          break;
        default:
          check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
          break;
        }
      break;

    case DNP1:
      switch (z_discretisation) {
        case LGP0:
          paire2D=DNP1xLGP0;
          break;
        case LGP1:
          paire2D=DNP1xLGP1;
          break;
        case LGP2:
          paire2D=DNP1xLGP2;
          break;
        default:
          check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
          break;
        }
      break;

    case CQP1:
      switch (z_discretisation) {
        case CQP0:
          paire2D=CQP1xCQP0;
          break;
        default:
          check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
          break;
       }
      break;

    case CQN1:
      switch (z_discretisation) {
        case CQP0:
          paire2D=CQN1xCQP0;
          break;
        default:
          check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
          break;
       }
      break;

    default:
      check_error(-1, "illicit u discretisation", __LINE__, __FILE__, 1);
      break;
   }
  return(paire2D);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  toponyms_t paire_map(void)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  toponyms_t Pmap;

  Pmap["LGP1xLGP0"]=LGP1xLGP0;
  Pmap["LGP0xLGP1"]=LGP0xLGP1;
  Pmap["LGP1xLGP1"]=LGP1xLGP1;
  Pmap["LGP1xLGP2"]=LGP1xLGP2;
  Pmap["NCP1xLGP0"]=NCP1xLGP0;
  Pmap["NCP1xLGP1"]=NCP1xLGP1;
  Pmap["NCP1xLGP2"]=NCP1xLGP2;
  Pmap["DGP1xLGP2"]=DGP1xLGP2;
  Pmap["DNP1xLGP0"]=DNP1xLGP0;
  Pmap["DNP1xLGP1"]=DNP1xLGP1;
  Pmap["DNP1xLGP2"]=DNP1xLGP2;

  return(Pmap);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int paire_dimension(mesh_t mesh,int paire, int *zdim,int *udim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *zdim=-1;
  *udim=-1;

  switch(paire) {
    case NCP1xLGP0:
/* *----------------------------------------------------------------------
      NCP1xLGP0 structure initialisation*/
      *zdim=mesh.ntriangles;
      *udim=mesh.nedges;
      break;

    case NCP1xLGP1:
/* *----------------------------------------------------------------------
      NCP1xLGP1 structure initialisation*/
      *zdim=mesh.nvtxs;
      *udim=mesh.nedges;
      break;

    case NCP1xLGP2:
/* *----------------------------------------------------------------------
      NCP1xLGP1 structure initialisation*/
      *zdim=mesh.nvtxs+mesh.nedges;
      *udim=mesh.nedges;
      break;

    case LGP0xLGP1:
/* *----------------------------------------------------------------------
      LGP1xLGP0 structure initialisation*/
      *zdim=mesh.nvtxs;
      *udim=mesh.ntriangles;
      break;

    case LGP1xLGP0:
/* *----------------------------------------------------------------------
      LGP1xLGP0 structure initialisation*/
      *zdim=mesh.ntriangles;
      *udim=mesh.nvtxs;
      break;

    case LGP1xLGP1:
/* *----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *zdim=mesh.nvtxs;
      *udim=mesh.nvtxs;
      break;

    case LGP1xLGP2:
/* *----------------------------------------------------------------------
      LGP1xLGP2 structure initialisation*/
      *zdim=mesh.nvtxs+mesh.nedges;
      *udim=mesh.nvtxs;
      break;

    case DGP1xLGP2:
/* *----------------------------------------------------------------------
      DGP1xLGP2 structure initialisation*/
      if(mesh.LGP2descriptor.NIbE==0) return(-1);
      *zdim=mesh.LGP2descriptor.nnodes;
      if(mesh.DGP1descriptor.NIbE==0) return(-1);
      *udim=mesh.DGP1descriptor.nnodes;
      break;

    case DNP1xLGP2:
/* *----------------------------------------------------------------------
      DNP1xLGP2 structure initialisation*/
      if(mesh.LGP2descriptor.NIbE==0) return(-1);
      *zdim=mesh.LGP2descriptor.nnodes;
      if(mesh.DNP1descriptor.NIbE==0) return(-1);
      *udim=mesh.DNP1descriptor.nnodes;
      break;

    default:
      check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
      break;
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int vector_dimension(mesh_t mesh,int discretisation, int *dim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  switch(discretisation) {
    case LGP0:
/* *----------------------------------------------------------------------
      LGP0 structure initialisation*/
      *dim=mesh.ntriangles;
      break;

    case LGP1:
/* *----------------------------------------------------------------------
      LGP1 structure initialisation*/
      *dim=mesh.nvtxs;
      break;

    case NCP1:
/* *----------------------------------------------------------------------
      NCP1 structure initialisation*/
      *dim=mesh.nedges;
      break;

    case DGP1:
/* *----------------------------------------------------------------------
      DGP1 structure initialisation*/
      if(mesh.DGP1descriptor.NIbE==0) return(-1);
      *dim=mesh.DGP1descriptor.nnodes;
      break;

    case LGP2:
/* *----------------------------------------------------------------------
      LGP2 structure initialisation*/
      if(mesh.LGP2descriptor.NIbE==0) return(-1);
      *dim=mesh.LGP2descriptor.nnodes;
      break;

    case DNP1:
/* *----------------------------------------------------------------------
      DNP1 structure initialisation*/
      if(mesh.DNP1descriptor.NIbE==0) return(-1);
      *dim=mesh.DNP1descriptor.nnodes;
      break;

    case CQP0:
/**----------------------------------------------------------------------
      */
      if(mesh.CQP0descriptor.NIbE==0) return(-1);
      *dim=mesh.CQP0descriptor.nnodes;
      break;

    case CQP1:
/**----------------------------------------------------------------------
      */
      if(mesh.CQP1descriptor.NIbE==0) return(-1);
      *dim=mesh.CQP1descriptor.nnodes;
      break;

    case CQN1:
/**----------------------------------------------------------------------
      */
      if(mesh.CQN1descriptor.NIbE==0) return(-1);
      *dim=mesh.CQN1descriptor.nnodes;
      break;

    default:
      check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
      break;
    }
  return(0);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int vector_dimension(mesh_t mesh, int discretisation, int *dim, int context)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   switch(context) {
//     case CONTEXT_LOCAL:
//       status=vector_dimension(mesh,discretisation,dim);
//       break;
//     case CONTEXT_GLOBAL:
//       if(mesh.origin==0) {
//         status=vector_dimension(mesh,discretisation,dim);
//         }
//       else {
// //        printf("paire_dimension, global context: %d\n",mesh.origin);
//         status=vector_dimension(*(mesh.origin),discretisation,dim);
//         }
//       break;
//     default:
//       check_error(-1, "unknown context", __LINE__, __FILE__, 1);
//       break;
//     }
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int solved_dimension(mesh_t mesh,int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  discretisation_t descriptor;

  descriptor= get_descriptor(mesh, discretisation);

#ifdef PARALLEL
/* *----------------------------------------------------------------------------
  compute nb nodes locally resolved */
  nsolved=0;
  for (i=0; i <descriptor.nnodes; i++) {
    if ( (descriptor.nodes[i].rank) == 0) {
      nsolved++;
      }
    }

  return(nsolved);
#else
  return(descriptor.nnodes);
#endif
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int paire_dimension(mesh_t mesh,int *zdim, int *udim)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//
//   status=paire_dimension(mesh,gSequentialPaire2D,zdim,udim);
// }

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int paire_dimension_contextual(mesh_t mesh,int *zdim, int *udim, int context)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   switch(context) {
//     case CONTEXT_LOCAL:
//       status=paire_dimension(mesh,gSequentialPaire2D,zdim,udim);
//       break;
//     case CONTEXT_GLOBAL:
//       if(mesh.origin==0) {
//         status=paire_dimension(mesh,gSequentialPaire2D,zdim,udim);
//         }
//       else {
// //        printf("paire_dimension, global context: %d\n",mesh.origin);
//         status=paire_dimension(*(mesh.origin),gSequentialPaire2D,zdim,udim);
//         }
//       break;
//     default:
//       check_error(-1, "unknown context", __LINE__, __FILE__, 1);
//       break;
//     }
// }

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int paire_discretisation_id(int paire, int *z_discretisation,int *u_discretisation,int *g_discretisation)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   status=paire_discretisation_id(paire, z_discretisation, u_discretisation);
//    if(g2D_discretisation==AUTO) {
//      *g_discretisation=gradient_natural_discretisation(u2D_discretisation);
//      }
//    if(g2D_discretisation==INTRINSIC) {
//      *g_discretisation=gradient_natural_discretisation(u2D_discretisation);
//      }
//    else {
//      *g_discretisation=g2D_discretisation;
//      }
// }

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int paire_dimension(mesh_t mesh,int *zdim,int *udim, int *gdim)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int g_discretisation, status;
//
//    status=paire_dimension(mesh,gSequentialPaire2D,zdim,udim);
//    if(g2D_discretisation==AUTO) {
//      g_discretisation=gradient_natural_discretisation(u2D_discretisation);
//      }
//    if(g2D_discretisation==INTRINSIC) {
//      g_discretisation=gradient_natural_discretisation(u2D_discretisation);
//      }
//    else {
//      g_discretisation=g2D_discretisation;
//      }
//    status=vector_dimension(mesh,g_discretisation,gdim);
// }

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int paire_order(int *zorder, int *uorder)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
//
//   Return dynamics equation order
//
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
// {
//   int status;
//   int integration;
//   int j, n;
//
// /* *----------------------------------------------------------------------
//   set matrix dynamical coefficients*/
//   switch (we_type) {
//     case GWE_CLX:
// /* *----------------------------------------------------------------------
//       Generalized wave equation, Lynch & Gray [1979] formulation*/
//       *zorder=1;
//       *uorder=1;
//       break;
//
//     case GWE_XXX:
// /* *----------------------------------------------------------------------
//       Generalized wave equation, revisited Lynch & Gray [1979] formulation*/
//       *zorder=2;
//       *uorder=1;
//       break;
//
// /* *----------------------------------------------------------------------
//     Explicit wave equation*/
//     case WE_EXPLICIT:
//       *zorder=1;
//       *uorder=1;
//       break;
//
// /* *----------------------------------------------------------------------
//     Implicit wave equation*/
//     case WE_IMPLICIT:
//       *zorder=1;
//       *uorder=1;
//       break;
//
//     default:
//       check_error(-1, "illegal mode", __LINE__, __FILE__, 1);
//       break;
//   }
//
//   return (0);
// }


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int triplet_discretisation_id(int triplet, int *w_discretisation, int *u_discretisation, int *t_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  switch(triplet) {
    case LGP0xQLP1xLGP1:
      *w_discretisation=QLP1;
      *u_discretisation=LGP0;
      *t_discretisation=LGP1;
      break;
    case NCP1xQLP0xLGP0:
      *w_discretisation=QLP0;
      *u_discretisation=NCP1;
      *t_discretisation=LGP0;
      break;

    default:
      check_error(-1, "unknown triplet (not implemented?)", __LINE__, __FILE__, 1);
      break;
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int triplet_dimension(mesh_t mesh, int triplet, int *w_dim, int *u_dim, int *t_dim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  int w_discretisation, u_discretisation, t_discretisation;
  
  status=triplet_discretisation_id(triplet, &w_discretisation, &u_discretisation, &t_discretisation);
  
  status= vector_dimension(mesh, w_discretisation, w_dim);
  status= vector_dimension(mesh, u_discretisation, u_dim);
  status= vector_dimension(mesh, t_discretisation, t_dim);
 
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int triplet_definition(int u_discretisation,int w_discretisation,int t_discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int triplet;

  switch (u_discretisation) {
    case LGP0:
      switch (w_discretisation) {
//         case QLP0:
//           switch (t_discretisation) {
//             case LGP0:
//                triplet = LGP0xQLP0xLGP0;
//                break;
//             case LGP1:
//                triplet = LGP0xQLP0xLGP1;
//                break;
//             case QLP1:
//                triplet = LGP0xQLP0xQLP1;
//                break;
//             }
//           break;
        case QLP1:
          switch (t_discretisation) {
//             case LGP0:
//                triplet = LGP0xQLP1xLGP0;
//                break;
            case LGP1:
               triplet = LGP0xQLP1xLGP1;
               break;
//             case QLP1:
//                triplet = LGP0xQLP1xQLP1;
//                break;
            }
          break;
        default:
          check_error(-1, "discretisation not (yet?) implemented", __LINE__, __FILE__, 1);
          break;
        }
      break;
    case NCP1:
      switch (w_discretisation) {
        case QLP0:
          switch (t_discretisation) {
            case LGP0:
               triplet = NCP1xQLP0xLGP0;
               break;
            case LGP1:
               triplet = NCP1xQLP0xLGP1;
               break;
            case QLP1:
               triplet = NCP1xQLP0xQLP1;
               break;
            }
          break;
        case QLP1:
          switch (t_discretisation) {
            case LGP0:
               triplet = NCP1xQLP1xLGP0;
               break;
            case LGP1:
               triplet = NCP1xQLP1xLGP1;
               break;
            case QLP1:
               triplet = NCP1xQLP1xQLP1;
               break;
            }
          break;
        default:
          check_error(-1, "discretisation not (yet?) implemented", __LINE__, __FILE__, 1);
          break;
        }
      break;
    default:
      check_error(-1, "discretisation not (yet?) implemented", __LINE__, __FILE__, 1);
      break;
    }

  return(triplet);
}


/**
    return discretisation paire id from its name
    \date 18 May 2011
    \author Laurent Roblou
 */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int paire_discretisation_init(string discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  if(discretisation.compare("NCP1xLGP0") == 0) {
    return(NCP1xLGP0);
  }
  if(discretisation.compare("NCP1xLGP1") == 0) {
    return(NCP1xLGP1);
  }
  if(discretisation.compare("NCP1xLGP2") == 0) {
    return(NCP1xLGP2);
  }
  if(discretisation.compare("LGP1xLGP0") == 0) {
    return(LGP1xLGP0);
  }
  if(discretisation.compare("NCP1xLGP0") == 0) {
    return(NCP1xLGP0);
  }
  if(discretisation.compare("LGP1xLGP1") == 0) {
    return(LGP1xLGP1);
  }
  if(discretisation.compare("LGP1xLGP2") == 0) {
    return(LGP1xLGP2);
  }
  if(discretisation.compare("DGP1xLGP0") == 0) {
    return(DGP1xLGP0);
  }
  if(discretisation.compare("DGP1xLGP0") == 0) {
    return(DGP1xLGP0);
  }
  if(discretisation.compare("DGP1xLGP1") == 0) {
    return(DGP1xLGP1);
  }
  if(discretisation.compare("DGP1xLGP2") == 0) {
    return(DGP1xLGP2);
  }
  if(discretisation.compare("DNP1xLGP0") == 0) {
  return(DNP1xLGP0);
}
  if(discretisation.compare("DNP1xLGP1") == 0) {
  return(DNP1xLGP1);
}
  if(discretisation.compare("DNP1xLGP2") == 0) {
  return(DNP1xLGP2);
}
  if(discretisation.compare("LGP0xLGP1") == 0) {
    return(LGP0xLGP1);
  }
 
  return(-1);

}


/**
    return discretisation single ids from a paire
    \date 19 May 2011
    \author F. Lyard
 */

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
//
//   int paire_discretisation_id(int paire, int *z_discretisation, int *u_discretisation)
//
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//
//   switch(paire) {
//     case NCP1xLGP0:
// /* *----------------------------------------------------------------------
//       NCP1xLGP0 structure initialisation*/
//       *z_discretisation=LGP0;
//       *u_discretisation=NCP1;
//       break;
//
//     case NCP1xLGP1:
// /* *----------------------------------------------------------------------
//       NCP1xLGP1 structure initialisation*/
//       *z_discretisation=LGP1;
//       *u_discretisation=NCP1;
//       break;
//
//     case NCP1xLGP2:
// /* *----------------------------------------------------------------------
//       NCP1xLGP1 structure initialisation*/
//       *z_discretisation=LGP2;
//       *u_discretisation=NCP1;
//       break;
//
//     case LGP0xLGP1:
// /* *----------------------------------------------------------------------
//       LGP1xLGP0 structure initialisation*/
//       *z_discretisation=LGP1;
//       *u_discretisation=LGP0;
//       break;
//
//     case LGP1xLGP0:
// /* *----------------------------------------------------------------------
//       LGP1xLGP0 structure initialisation*/
//       *z_discretisation=LGP0;
//       *u_discretisation=LGP1;
//       break;
//
//     case LGP1xLGP1:
// /* *----------------------------------------------------------------------
//       LGP1xLGP1 structure initialisation*/
//       *z_discretisation=LGP1;
//       *u_discretisation=LGP1;
//       break;
//
//     case LGP1xLGP2:
// /* *----------------------------------------------------------------------
//       LGP1xLGP2 structure initialisation*/
//       *z_discretisation=LGP2;
//       *u_discretisation=LGP1;
//       break;
//
//     case DGP1xLGP2:
// /* *----------------------------------------------------------------------
//       LGP1xLGP1 structure initialisation*/
//       *z_discretisation=LGP2;
//       *u_discretisation=DGP1;
//       break;
//
// //     case DNP1xLGP2:
// // /* *----------------------------------------------------------------------
// //       LGP1xLGP1 structure initialisation*/
// //       *z_discretisation=LGP2;
// //       *u_discretisation=DNP1;
// //       break;
//
//     default:
//       check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
//       break;
//   }
//   return(0);
// }



/**
    return discretisation single dimensions from a paire
    \date 19 May 2011
    \author F. Lyard
 */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int paire_dimension(const mesh_t& mesh, int paire, int *zdim,int *udim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  *zdim=-1;
  *udim=-1;

  switch(paire) {
    case NCP1xLGP0:
/**----------------------------------------------------------------------
      NCP1xLGP0 structure initialisation*/
      *zdim=mesh.ntriangles;
      *udim=mesh.nedges;
      break;

    case NCP1xLGP1:
/**----------------------------------------------------------------------
      NCP1xLGP1 structure initialisation*/
      *zdim=mesh.nvtxs;
      *udim=mesh.nedges;
      break;

    case NCP1xLGP2:
/**----------------------------------------------------------------------
      NCP1xLGP1 structure initialisation*/
      *zdim=mesh.nvtxs+mesh.nedges;
      *udim=mesh.nedges;
      break;

    case LGP0xLGP1:
/**----------------------------------------------------------------------
      LGP1xLGP0 structure initialisation*/
      *zdim=mesh.nvtxs;
      *udim=mesh.ntriangles;
      break;

    case LGP1xLGP0:
/**----------------------------------------------------------------------
      LGP1xLGP0 structure initialisation*/
      *zdim=mesh.ntriangles;
      *udim=mesh.nvtxs;
      break;

    case LGP1xLGP1:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      *zdim=mesh.nvtxs;
      *udim=mesh.nvtxs;
      break;

    case LGP1xLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      if(mesh.LGP2descriptor.NIbE==0) return(-1);
      *zdim=mesh.LGP2descriptor.nnodes;
      *udim=mesh.nvtxs;
      break;

    case DGP1xLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      if(mesh.LGP2descriptor.NIbE==0) return(-1);
      *zdim=mesh.LGP2descriptor.nnodes;
      if(mesh.DGP1descriptor.NIbE==0) return(-1);
      *udim=mesh.DGP1descriptor.nnodes;
      break;

    case DNP1xLGP1:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      if(mesh.LGP1descriptor.NIbE==0) return(-1);
      *zdim=mesh.LGP1descriptor.nnodes;
      if(mesh.DNP1descriptor.NIbE==0) return(-1);
      *udim=mesh.DNP1descriptor.nnodes;
      break;

    case DNP1xLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      if(mesh.LGP2descriptor.NIbE==0) return(-1);
      *zdim=mesh.LGP2descriptor.nnodes;
      if(mesh.DNP1descriptor.NIbE==0) return(-1);
      *udim=mesh.DNP1descriptor.nnodes;
      break;

    case LGP2xLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      if(mesh.LGP2descriptor.NIbE==0) return(-1);
      *zdim=mesh.LGP2descriptor.nnodes;
      if(mesh.LGP2descriptor.NIbE==0) return(-1);
      *udim=mesh.LGP2descriptor.nnodes;
      break;

    case IPGxLGP2:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      if(mesh.LGP2descriptor.NIbE==0) return(-1);
      *zdim=mesh.LGP2descriptor.nnodes;
      if(mesh.IPGdescriptor.NIbE==0) return(-1);
      *udim=mesh.IPGdescriptor.nnodes;
      break;

    case CQP1xCQP0:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      if(mesh.CQP0descriptor.NIbE==0) return(-1);
      *zdim=mesh.CQP0descriptor.nnodes;
      if(mesh.CQP1descriptor.NIbE==0) return(-1);
      *udim=mesh.CQP1descriptor.nnodes;
      break;

    case CQN1xCQP0:
/**----------------------------------------------------------------------
      LGP1xLGP1 structure initialisation*/
      if(mesh.CQP0descriptor.NIbE==0) return(-1);
      *zdim=mesh.CQP0descriptor.nnodes;
      if(mesh.CQN1descriptor.NIbE==0) return(-1);
      *udim=mesh.CQN1descriptor.nnodes;
      break;

    default:
      check_error(-1, "unknown paire", __LINE__, __FILE__, 1);
      break;
    }
  return(0);
}


/**
    return discretisation dimensions from a paire
    \date 27 July 2011
    \author Laurent Roblou
 */
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

paire_t paire_discretisation_ids(string FEdescription)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int paire = paire_discretisation_init(FEdescription);
  
  int z_discretisation = 0, u_discretisation = 0;
  (void) paire_discretisation_id(paire, &z_discretisation, &u_discretisation);

  paire_t ids(z_discretisation, u_discretisation);
  return(ids);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

mesh_t& mesh_t::initialize_descriptors(string FePaire, int compute_massmatrix)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int paire = paire_discretisation_init(FePaire);
  int status;
      
  int z_discretisation = 0, u_discretisation = 0;
  (void) paire_discretisation_id(paire, &z_discretisation, &u_discretisation);
    
  status=discretisation_init(this, z_discretisation, compute_massmatrix); //!< \todo bugfix: Invalid write of size 4: initLGP2(mesh_t, int) (discretisation-initialise.cpp:512)
    
  status=discretisation_init(this, u_discretisation, compute_massmatrix);

  return(*this);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

mesh_t& mesh_t::initialize_descriptors(int FePaire, int compute_massmatrix)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int z_discretisation = 0, u_discretisation = 0;
  (void) paire_discretisation_id(FePaire, &z_discretisation, &u_discretisation);
    
  (void) discretisation_init(this, z_discretisation, compute_massmatrix); //!< \todo bugfix: Invalid write of size 4: initLGP2(mesh_t, int) (discretisation-initialise.cpp:512)
    
  (void) discretisation_init(this, u_discretisation, compute_massmatrix);

  return(*this);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_element2neighbours (mesh_t mesh, discretisation_t descriptor, int nnghbs_max)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  build the nodes' neighbours list from the elements table

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i,k,l,m,n,n1,n2;
  int **nghbs,acquire;
  
  nghbs=new int*[descriptor.nnodes];
  
  for(n = 0; n < descriptor.nnodes; n++) {
    descriptor.nodes[n].nnghbs=0;
    nghbs[n]=new int[nnghbs_max];
    for(i=0;i<nnghbs_max;i++) nghbs[n][i]=-1;
    }
    
/**----------------------------------------------------------------------------
  create neighbours connectivity (for each computational nodes)*/
  for(m = 0; m < mesh.ntriangles; m++) {
    for(k=0;k<descriptor.nnpe;k++) {
      n1=descriptor.NIbE[m][k];
      for(l=0;l<descriptor.nnpe;l++) {
        if(k==l) continue;
        n2=descriptor.NIbE[m][l];
        acquire=1;
        for(i=0;i<descriptor.nodes[n1].nnghbs;i++) {
          if(nghbs[n1][i]==n2) {
            acquire=0;
            break;
            }
          }
        if(acquire==1) {
          nghbs[n1][i]=n2;
          descriptor.nodes[n1].nnghbs++;
          if(descriptor.nodes[n1].nnghbs>=nnghbs_max) {
            printf("node %d (%d in element %d) : \n",n1,k,m);
            fe_print(mesh, m);
            check_error(-1, "discretisation, neigbours overflow", __LINE__, __FILE__, 1);
            }
          }
        }
      }
    }
  for(n = 0; n < descriptor.nnodes; n++) {
    descriptor.nodes[n].nghbs=new int[descriptor.nodes[n].nnghbs];
    for(i=0;i<descriptor.nodes[n].nnghbs;i++) descriptor.nodes[n].nghbs[i]=nghbs[n][i];
    }

  delete[] nghbs;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_element2neighbours (mesh_t mesh, discretisation_t descriptor)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  build the nodes' neighbours list from the elements table

@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int status;
/**----------------------------------------------------------------------------
  max number of neighbours predictor, discretisation-dependent...*/
  int nnghbs_max;
  
  
  switch(descriptor.type) {
    case LGP1:
      nnghbs_max=mesh.nnghm;
      status=fe_element2neighbours (mesh, descriptor,nnghbs_max);
      break;
    case DGP1:
      nnghbs_max=3;
      status=fe_element2neighbours (mesh, descriptor,nnghbs_max);
      break;
    case NCP1:
      nnghbs_max=4;
      status=fe_element2neighbours (mesh, descriptor,nnghbs_max);
      break;
    case DNP1:
      nnghbs_max=3;
      status=fe_element2neighbours (mesh, descriptor,nnghbs_max);
      break;
    case LGP2:
      nnghbs_max=mesh.nnghm*3+2;
      status=fe_element2neighbours (mesh, descriptor,nnghbs_max);
      break;
    case DGP2:
      nnghbs_max=6;
      status=fe_element2neighbours (mesh, descriptor,nnghbs_max);
      break;
    default:
      status=-1;
      break;
  }

  return(0);
}


