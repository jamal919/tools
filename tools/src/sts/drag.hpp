
#ifndef DRAG_HPP
#define DRAG_HPP

#include "exceptions.hpp"

void karman_friction(state2d_t, parameter_t, int) ;
void karman_Cd(parameter_t, double *, int, float *) ;
void karman_Cd3D(parameter_t, discretisation_t, int, float *) ;
//void karman_frictionLG(double, state2d_t, parameter_t, double, int) ;

#endif
