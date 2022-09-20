
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief discretisation classes and function declarations
*/
/*----------------------------------------------------------------------------*/

#ifndef DISCRETISATION_H
#define DISCRETISATION_H

#define PerColumns 0
#define PerLayers  1

#define Order3D PerColumns

#define  OFFSET 1

#include "matrix.h" /* for hypermatrix_t and distribution_t */


class node_t {
private :
public :
  double lon, lat, h;           /* longitude, latitude, mean depth    */
  double *zlevels, *sigma;       /* depth levels,sigma levels          */
  double c, s, weight;          /* cos/sin latitude                   */
  double g;                     /* gravitational constant             */
  int nlevels;                  /* sigma levels                       */
  int *nghbs,nnghbs;            /* neighbours (vertices) list         */
  int *vtces,nvtces;            /* neighbours (vertices) list         */
  int *edges,nedges;            /* connected edges list               */
  int *elmts,nelmts;            /* elements list                      */
  int code;                     /* boundary code for nodee            */
  int child,parent,ancestro;    /* index in imbricated mesh           */
  int globalized;               /* index in partitioned mesh          */
  int rank;                     /*                                    */

  node_t() {
    lon=9999.999;lat=9999.999;
    c=9999.999,s=9999.999;
    weight=0;
    zlevels=sigma=NULL;
    nghbs=0;nnghbs=-1;
    vtces=0;nvtces=-1;
    edges=0;nedges=-1;
    elmts=0;nelmts=-1;
    code=-1;
    child=-1;parent=-1;ancestro=-1;
    rank=-1;
    }

  void destroy() {
    if(nghbs!=0) {
      delete[] nghbs;
      nghbs=0;
      };
    nnghbs=-1;
    if(vtces!=0) {
      delete[] vtces;
      vtces=0;
      };
    nvtces=-1;
    if(edges!=0) {
      delete[] edges;
      edges=0;
      };
    nedges=-1;
    if(elmts!=0) {
      delete[] elmts;
      elmts=0;
      };
    nelmts=-1;
    if(sigma!=0) {
      delete[] sigma;
      sigma=0;
      };
    if(zlevels!=0) {
      delete[] zlevels;
      zlevels=0;
      };
    }
  };

class element__t {
private :
public :
  int *node;            /* neighbours (vertices) list         */
  element__t() {
    node=0;
    }
  };

  
class discretisation_t {
private :
public :
  int nelmts;                   /* number of elements                   */
  int nnodes;                   /* number of computational node         */
  int type;                     /* discretisation type                  */
  int nnpe;                     /* number of nodes per element          */
  int **NIbE;                   /* node index of element                */
  element__t *elements;         /*                                      */
  hypermatrix_t massmatrix;     /* mass matrix                          */
  node_t *nodes;                /* node structure                       */
  distribution_t distributor;   /* node distribution information (MPI)  */
  /* STS: NO periodic */

  discretisation_t() {
    nelmts=0;
    nnodes=0;
    type=-1;
    nnpe=0;
    NIbE=0;
    elements=0;
    nodes=0;
    }

  void destroy() {
    if(NIbE!=0) {
      int n;
      for(n=0;n<nelmts;n++) {
        delete[] NIbE[n];
        }
      delete[] NIbE;
      NIbE=0;
      };
    if(nodes!=0) {
      int n;
      for(n=0;n<nnodes;n++) {
        nodes[n].destroy();
        }
      delete[] nodes;
      nodes=0;
      };
    massmatrix.destroy();
    nnodes=0;
    type=0;
    nnpe=0;
    /* STS: NO periodic */
    }

};

// extern int structured3DIndex(int indexH, int indexV, int hdim, int vdim);
// extern int structured3DIndex(int indexH, int indexV, int indexB, int hdim, int vdim, int bdim);
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static inline int structured3DIndex(int indexH, int indexV, int hdim, int vdim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------

  return an implicit 3D node number as a vertical extension of a 2D (horizontal)
  numbering.
  
  indexH : horizontal node number
  indexV : layer
  
  used for monovariate nodes (typically p or tracer nodes)
 
-----------------------------------------------------------------------------*/
{
  switch (Order3D) {
    case PerColumns:
      return(indexH*vdim+indexV);
    case PerLayers:
      return(indexV*hdim+indexH);
    default:
      TRAP_ERR_EXIT(-1, "unknown 3D per levels/layers ordering\n");
   }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

static inline int structured3DIndex(int indexH, int indexV, int indexB, int hdim, int vdim, int bdim)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/**----------------------------------------------------------------------------

  return an implicit 3D node number as a vertical extension of a 2D (horizontal)
  numbering.
  
  indexH : horizontal node number
  indexV : layer
  
  used for multivariate nodes (typically [u,v] nodes)
  
-----------------------------------------------------------------------------*/
  switch (Order3D) {
    case PerColumns:
      return(indexB*(hdim*vdim)+indexH*vdim+indexV);
    case PerLayers:
      return(indexB*(hdim*vdim)+indexV*hdim+indexH);
    default:
      TRAP_ERR_EXIT(-1, "unknown 3D per levels/layers ordering\n");
   }
  return(0);
}

extern int gradient_natural_discretisation(int discretisation);

//#define indexLGP0xLGP1_2D3D(pointer,K,L,n,idx,k,l) (pointer[n]*K*L+k*L*3+idx*L+l)
#define indexLGP0xLGP1_2D3D(pointer,K,L,n,idx,k,l) ((pointer[n]+idx)*K*L+k*L+l)

#endif /* DISCRETISATION_H */
