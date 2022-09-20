#ifndef DRY_H
#   define DRY_H

#   ifdef MAIN_SOURCE
#      define PREFIX
#   else
#      define PREFIX extern
#   endif

PREFIX double shallowest;
PREFIX int    *gzNodeWet, *gTriangleWet;

//#   ifdef WEDBLE
PREFIX double *ActualMatrix[NLVL];
PREFIX double *DryingMatrix[NLVL];
// #   else
// PREFIX float  *ActualMatrix[NLVL];
// PREFIX float  *DryingMatrix[NLVL];
// #   endif

#endif /* DRY_H */
