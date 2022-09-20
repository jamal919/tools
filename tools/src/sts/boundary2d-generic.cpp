
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

#include "tugo-prototypes.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int VelocityBCs_init(specification_obsolete_t spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i;
  double size;

  if (tugo_cfg->boundaries->SolidCondition.face_value() == "FREE-SLIP") {
    size=0.;
    }
  else if (tugo_cfg->boundaries->SolidCondition.face_value() == "CONDITIONAL-NO-SLIP") {
    size=tugo_cfg->boundaries->NoSlipSize.numerical_value();
    }
  
  for(i = 0; i < spec.nnodes; i++) {
    spec.fluxes[i].option=FREE_SLIP;
    if(spec.data[i].L < size) {
      spec.fluxes[i].option=NO_SLIP;
      }
    }

  return(0);
}
