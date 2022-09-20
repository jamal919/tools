
/**************************************************************************

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

***************************************************************************/

#include "tugo-prototypes.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void update_U_elevation(mesh_t & mesh, state2d_t & state2D,parameter_t & data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Compute H elevation at U nodes

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i, j, m, n, m1,m2;
  int status;

  status=fe_projection(mesh, state2D.H, z2D_discretisation, state2D.H__,u2D_discretisation);
  status=fe_projection(mesh, state2D.z, z2D_discretisation, state2D.z__,u2D_discretisation);

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void update_uNodes_elevation(mesh_t & mesh, state2d_t & state2D,parameter_t & data)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Compute H elevation at U nodes

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check :

  Note:


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
{
  int i, j, m, n, n1, n2;
  int status;


  switch (u2D_discretisation) {

    case LGP0:
/*----------------------------------------------------------------------
      DGP1 elevation discretisation */
      update_U_elevation( mesh, state2D, data);
      break;

      /* STS: NO LGP1 */

    case DGP1:
/*----------------------------------------------------------------------
      DGP1 elevation discretisation */
      update_U_elevation( mesh, state2D, data);
      break;

      /* STS: NO NCP1 */

    case DNP1:
/*----------------------------------------------------------------------
      LGP1 elevation discretisation */
      update_U_elevation( mesh,  state2D, data);
      break;

    case CQP0:
/*----------------------------------------------------------------------
      LGP1 elevation discretisation */
      update_U_elevation( mesh,  state2D, data);
      break;

    case CQP1:
/*----------------------------------------------------------------------
      LGP1 elevation discretisation */
      update_U_elevation( mesh,  state2D, data);
      break;

    case CQN1:
/*----------------------------------------------------------------------
      LGP1 elevation discretisation */
      update_U_elevation( mesh,  state2D, data);
      break;

    default:
      check_error(-1, "illegal elevation discretisation", __LINE__, __FILE__, 1);
      break;
    }
}
