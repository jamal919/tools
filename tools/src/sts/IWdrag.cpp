
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
#include "constants.h"
#include "keywords.hpp"
#include "map.h"
#include "grd.h"
#include "geo.h"
#include "maths.h"

#include "polygons.hpp"

extern float *tide_prediction_constituents(double, spectrum_t, hconstant_t, int);
extern float tide_prediction(double, spectrum_t, hconstant_t, int);

extern  void ellipse_parameter02(float Au, float Gu, float Av, float Gv,float *a,float *b,float *c, float *e,float *dir);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void initialize_IWdrag(mesh_t & mesh, parameter_t & data, int nndes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  using namespace Keywords;
  int k,m,n,status;
  int count;
  char filename[1024], rootname[1024];
  float mask;
  extern int load_EigenModes(const char *filename, parameter_t & data, int nndes, float *celerity, float *N, float *spec);

  fprintf(echo, "xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Internal wave drag section\n");

  wdc1 = tugo_cfg->dissipation->InternDragSlope.numerical_value<double>() ;
  fprintf(echo, " Wave drag coefficient (slope) : %f\n", wdc1);

  wdc2 = tugo_cfg->dissipation->InternDragRugosity.numerical_value<double>() ;
  fprintf(echo, " Wave drag coefficient (rugosity) : %f\n", wdc2);

  gwdHmin = tugo_cfg->dissipation->InternDragHmin.numerical_value<double>() ;
  fprintf(echo, " Wave drag Hmin :- %f\n", gwdHmin);

  gwdHmax = tugo_cfg->dissipation->InternDragHmax.numerical_value<double>() ;
  fprintf(echo, " Wave drag Hmax : %f\n", gwdHmax);

  wdAlgorithm = tugo_cfg->dissipation->InternDragAlgorithm.numerical_value<int>() ;
  fprintf(echo, " Wave drag algoritm  : %d\n", wdAlgorithm);
  
/**----------------------------------------------------------------------------
  to keep consistent with older versions */
  for (n=0; n < nndes; n++) data.celerity[n] = 3.3;
  
/**----------------------------------------------------------------------------
  set flag to allow for multiple entries */
  for (n=0; n < nndes; n++) data.wdc1[n]=-1;
    
  if (tugo_cfg->dissipation->IWD01Polygons.face_value() != (string) "NONE") {
/**----------------------------------------------------------------------------
    initializing by polygons */
    /* STS: NO space-variable setting */
    }

  if (tugo_cfg->dissipation->IWD01ValuesByRegions.face_value() != (string) "NONE") {
/**----------------------------------------------------------------------------
    initializing by regions */
    /* STS: NO space-variable setting */
    }
  
/**----------------------------------------------------------------------------
  initializing by uniform coefficient */
  sprintf(filename,"/home/softs/data/climatology/levitus/levitus94.vertical-modes.nc");
//     status=load_EigenModes(filename, data, nndes, data.celerity, data.Nbar, &mask); 
//     do {
//       count=fill(mesh, data.celerity, mask, nndes);
//       } while(count!=0);
//     do {
//       count=fill(mesh, data.Nbar, mask, nndes);
//       } while(count!=0);
  for (n=0; n < nndes; n++) 
    if(data.wdc1[n]==-1) data.wdc1[n] = wdc1; 
  
/**----------------------------------------------------------------------------
  to keep consistent with older versions */
  for (n=0; n < nndes; n++) data.wdc1[n] /= 100.;

  wave_drag = ((wdc1 != 0.) || (wdc2 != 0.) || (tugo_cfg->dissipation->IWD01Polygons.face_value() != "NONE"));
}
