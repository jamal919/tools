

/**************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

#define MAIN_SOURCE

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tools-structures.h"

#include "geo.h"
#include "fe.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "functions.h"
#include "matrix.h"
#include "solverlib.h"

#include "statistic.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Mesh2DElasticDistorsionT(mesh_t & mesh, double* K, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
 * 
 * aimed to provide a method for smooth mesh distorsion by computing an elasticity
 * matrix from the original mesh that can be used to recompute vertices positions
 * in the distorded mesh after boundary edges modification
 * 
 * here it is coded for triangle meshes, easily adaptable to any other types of
 * meshes (quadrangles,...)
 * 
 * Next step will be a displacement operator to fit mesh boundaries with a given
 * line or polygon.
 * 
 * Boundary line-wise elasticity to be implemented:
 * 
 *   - 1 fixed vertex per boundary
 *   - equilibrium condition for the others: k1 * L1 + k2 * L2 =0
 * 
 *------------------------------------------------------------------------------ 
  
  Compute edge geometrical elasticity
  -----------------------------------
  
  Constrained problem:
  
    n interior vertices equations (to be strictly implemented)

    min potential energy :  E=(sum on neighbours) ki(xn-xi)² -> (sum on neighbours) ki xn = (sum) ki xi
  
    min departure from 1 :  E=(sum on edges) (ki-1)²
  
  
  Unconstrained problem (Lagrange mutipliers):
    
    min E=(sum on edges) (ki-1)² + (sum on vertices) ln (sum on neighbours j) kj(xn-xj), where ki and ln are unknowns
  
    dE/dki -> 2(ki-1) + (sum) ln (xn-xi) (2 vertices n involved)
  
    dE/dli -> (sum on neighbours j) kj(x-xj) = 0 (nngh edges involved)
  
    boundary conditions:
  
    ki = 1 for boundary edges
  
    li = 0 for boundary vertices
  
  
  Reference : https://fr.wikipedia.org/wiki/Multiplicateur_de_Lagrange
  
------------------------------------------------------------------------------*/
{
  int e,k,l,m,n,n1,n2,status;
  int row,col;
  size_t pos;
  hypermatrix_t R, P;
  double *rhs, *weight,Energy,*energyP1;
  double *x,*y;
  int count;
  int nnz;
  bool force_duplicate=false;
  int verbose=0;
//   bool debug=false;
  
  size_t neq=2*mesh.nvtxs+mesh.nedges;
  
  int *vertex=new int[mesh.nvtxs],*edge=new int[mesh.nedges];
  
  double mask=-1;
  
  if(K==0) K=new double[mesh.nedges];
  
  status=discretisation_init(&mesh,LGP1,0);
  
  x=new double[mesh.nvtxs];
  y=new double[mesh.nvtxs];
  
  weight=new double[mesh.nedges];
  for(n=0;n<mesh.nedges; n++) weight[n]=1.0;

  for(n=0;n<mesh.nvtxs; n++) {
    x[n]=mesh.vertices[n].lon;
    y[n]=mesh.vertices[n].lat;
    }
  
  for(n=0;n<mesh.nvtxs; n++) vertex[n]=-1;
  
  neq=mesh.nvtxs;
  
  P.ordering=new ordering_t;
  
  P.ordering->nrows=neq;
  
  P.ordering->cardinal=new int[neq];
  P.ordering->pointer =new int[neq+1];
  
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
      P.ordering->cardinal[n]  =1;
      }
    else {
      P.ordering->cardinal[n]  =1+mesh.vertices[n].nngh;
      }
    }
  
  P.ordering->pointer[0]=0;
  for(n=0;n<neq; n++) {
    P.ordering->pointer[n+1]=P.ordering->pointer[n]+P.ordering->cardinal[n];
    }
  
  nnz=P.ordering->pointer[neq];
  P.ordering->incidence=new int[nnz];
  
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
      pos=P.ordering->pointer[n];
      P.ordering->incidence[pos]=n;
      }
    else {
      int count=0;
      pos=P.ordering->pointer[n];
      P.ordering->incidence[pos+count]=n;
      count++;
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        m=mesh.vertices[n].ngh[k];
        P.ordering->incidence[pos+count]=m;
        count++;
        }
      }
    }
    
  status=matrix_reorder(P.ordering);
  
  P.allocate();
  rhs=new double[neq];
  
  P.assign(0.0);
    
//   for(n=0;n<mesh.nvtxs; n++) {
//     if(mesh.vertices[n].code!=0) {
//       row=n;
//       pos=P.ordering->pos(row,row);
//       P.packed[pos] =1.0;
//       rhs[row]=x[n];
//       }
//     else {
//       for(k=0;k<mesh.vertices[n].nngh;k++) {
//         int m=mesh.vertices[n].ngh[k];
//         int e=fe_isedge(mesh,n,m);
//         double dx=x[n]-x[m];
//         row=n;
//         pos=P.ordering->pos(row,row);
//         P.packed[pos] += K[e];
// 	col=m;
//         pos=P.ordering->pos(row,col);
//         P.packed[pos] = -K[e];
//         rhs[row]=0.0;
//         }
//       }
//     }

/*------------------------------------------------------------------------------
  set position matrix coefficients */
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
      row=n;
      pos=P.ordering->pos(row,row);
      P.packed[pos] =1.0;
      }
    else {
      for(k=0;k<mesh.vertices[n].nngh;k++) {
        int m=mesh.vertices[n].ngh[k];
        int e=fe_isedge(mesh,n,m);
        double dx=x[n]-x[m];
        row=n;
        pos=P.ordering->pos(row,row);
        P.packed[pos] += K[e];
        col=m;
        pos=P.ordering->pos(row,col);
        P.packed[pos] = -K[e];
        }
      }
    }
      
  P.distributor=new distribution_t;
  P.distributor->context=SEQUENTIAL_COMPUTING;
      
  status=LinearSystem_initialise(&P, 0, 0, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
   
/*------------------------------------------------------------------------------
  solve for longitudes */
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
      row=n;
      rhs[row]=x[n];
      }
    else {
      row=n;
      rhs[row]=0.0;
      }
    }
  
  status=LinearSystem_solve(P, rhs);
  
  count=0;
  for(n=0;n<mesh.nvtxs; n++) {
    mesh.vertices[n].lon=rhs[n];
    double d=fabs(x[n]-rhs[n]);
    if(d>1.e-06) {
//       printf("edge %d: K=%lf\n",K[n]);
      count++;
      }
    }

/*------------------------------------------------------------------------------
  solve for latitudes */
  for(n=0;n<mesh.nvtxs; n++) {
    if(mesh.vertices[n].code!=0) {
      row=n;
      rhs[row]=y[n];
      }
    else {
      row=n;
      rhs[row]=0;
      }
    }
      
  status=LinearSystem_solve(P, rhs);
  
  count=0;
  for(n=0;n<mesh.nvtxs; n++) {
    mesh.vertices[n].lat=rhs[n];
    double d=fabs(y[n]-rhs[n]);
    if(d>1.e-06) {
//       printf("edge %d: K=%lf\n",K[n]);
      count++;
      }
    }

  P.destroy();
  
  delete[] rhs;
  delete[] weight;
  
  delete[] x;
  delete[] y;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Mesh1DElasticDistorsionT(mesh_t & mesh, const vector<plg_t> & polygons, const vector<plg_t> & target, vector<plg_t> & optimal, double* K, double Kmove, double ratio, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int e,k,l,m,n,n1,n2,p,status;
  int row,col;
  size_t pos;
  hypermatrix_t P;
  double *rhs;
  double *x1,*y1,*x2,*y2;
  int count;
  size_t neq, nvtxs, nedges;
  double K1,K2;
  bool force_duplicate=false;
  int verbose=0;
  
// #define PERIODIC
  
  K1=(1.-ratio)*Kmove;
  K2=ratio*Kmove;

  if(polygons.size()!=target.size()) return(-1);
  
  neq=0;
  for(p=0;p<polygons.size();p++) {
    if(polygons[p].npt!=target[p].npt) return(-1);
    neq+=polygons[p].npt;
    }
  
  for(p=0;p<polygons.size();p++) {
    plg_t tmp;
    optimal.push_back(tmp);
    optimal[p].duplicate(target[p]);
    }
  
  nvtxs=neq;
  
  int *vertex=new int[nvtxs];
      
  x1=new double[nvtxs];
  y1=new double[nvtxs];
  
  x2=new double[nvtxs];
  y2=new double[nvtxs];
  
  n=0;
  for(p=0;p<polygons.size();p++) {
    for(k=0;k<polygons[p].npt;k++) {
      x1[n]=polygons[p].t[k];
      y1[n]=polygons[p].p[k];
      x2[n]=target[p].t[k];
      y2[n]=target[p].p[k];
      n++;
      }
    }
  
  for(n=0;n<nvtxs; n++) vertex[n]=-1;
  
  neq=nvtxs;
  
  P.ordering=new ordering_t;
  
  P.ordering->nrows=neq;
  
  P.ordering->cardinal=new int[neq];
  P.ordering->pointer =new int[neq+1];
  
  n=0;
  for(p=0;p<polygons.size();p++) {
#ifndef PERIODIC
    P.ordering->cardinal[n] =1;
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
      P.ordering->cardinal[n] =3;
      n++;
      }
    P.ordering->cardinal[n] =1;
    n++;
#else
    P.ordering->cardinal[n] =2;
    n++;
    for(k=1;k<polygons[p].npt;k++) {
      P.ordering->cardinal[n] =3;
      n++;
      }
#endif
    }
      
  P.ordering->pointer[0]=0;
  for(n=0;n<neq; n++) {
    P.ordering->pointer[n+1]=P.ordering->pointer[n]+P.ordering->cardinal[n];
    }
  
  size_t nnz=P.ordering->pointer[neq];
  P.ordering->incidence=new int[nnz];
  
  n=0;
  for(p=0;p<polygons.size();p++) {
#ifndef PERIODIC
    pos=P.ordering->pointer[n];
    P.ordering->incidence[pos]=n;
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
      pos=P.ordering->pointer[n];
      P.ordering->incidence[pos]  =n-1;
      P.ordering->incidence[pos+1]=n;
      P.ordering->incidence[pos+2]=n+1;
      n++;
      }
    pos=P.ordering->pointer[n];
    P.ordering->incidence[pos]=n;
    n++;
#else
    int first=n;
    pos=P.ordering->pointer[n];
    P.ordering->incidence[pos]=n;
    P.ordering->incidence[pos+1]=n+polygons[p].npt-1;
    n++;
    for(k=1;k<polygons[p].npt;k++) {
      pos=P.ordering->pointer[n];
      P.ordering->incidence[pos]  =n-1;
      P.ordering->incidence[pos+1]=n;
      int np1=first+(k+1) % polygons[p].npt;
      P.ordering->incidence[pos+2]=np1;
      n++;
      }
#endif
    }
              
  status=matrix_reorder(P.ordering);
   
  P.allocate();
  rhs=new double[neq];
  
  P.assign(0.0);
    
/*------------------------------------------------------------------------------
  set position matrix coefficients */
  n=0;
  for(p=0;p<polygons.size();p++) {
#ifndef PERIODIC
    row=n;
    pos=P.ordering->pos(row,row);
    P.packed[pos] =1.0;
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      row=n;
      e=n-1;
      pos=P.ordering->pos(row,row);
      P.packed[pos] += K[e];
      col=n-1;
      pos=P.ordering->pos(row,col);
      P.packed[pos] -= K[e];
      e=n;
      pos=P.ordering->pos(row,row);
      P.packed[pos] += K[e];
      col=n+1;
      pos=P.ordering->pos(row,col);
      P.packed[pos] -= K[e];
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      pos=P.ordering->pos(row,row);
      P.packed[pos] += K1+K2;
      n++;
      }
    row=n;
    pos=P.ordering->pos(row,row);
    P.packed[pos] =1.0;
    n++;
#else
    int first=n;
    row=n;
    pos=P.ordering->pos(row,row);
    P.packed[pos]   = 1.0;
    P.packed[pos+1] =-1.0;
    n++;
    for(k=1;k<polygons[p].npt;k++) {
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      row=n;
      e=n-1;
      pos=P.ordering->pos(row,row);
      P.packed[pos] += K[e];
      col=n-1;
      pos=P.ordering->pos(row,col);
      P.packed[pos] -= K[e];
      e=n;
      pos=P.ordering->pos(row,row);
      P.packed[pos] += K[e];
      int np1=first+(k+1) % polygons[p].npt;
      col=np1;
      pos=P.ordering->pos(row,col);
      P.packed[pos] -= K[e];
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      pos=P.ordering->pos(row,row);
      P.packed[pos] += K1+K2;
      P.packed[pos] += 2*K1;
      n++;
      }
#endif
    }
      
  P.distributor=new distribution_t;
  P.distributor->context=SEQUENTIAL_COMPUTING;
      
  status=LinearSystem_initialise(&P, 0, 0, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
      
/*------------------------------------------------------------------------------
  solve for longitudes */  
  n=0;
  for(p=0;p<polygons.size();p++) {
#ifndef PERIODIC
    row=n;
    rhs[row]=x2[n];
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      row=n;
      rhs[row]=K1*x1[n]+K2*x2[n];
      n++;
      }
    row=n;
    rhs[row]=x2[n];
    n++;
#else
    row=n;
    rhs[row]=0;
    n++;
    for(k=1;k<polygons[p].npt;k++) {
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      row=n;
      rhs[row]=K1*x1[n]+K2*x2[n];
      rhs[row]+=K1*x1[n-1]+K1*x1[n+1];
      n++;
      }
#endif
    }
  
  status=LinearSystem_solve(P, rhs);
  
  count=0;
  n=0;
  for(p=0;p<polygons.size();p++) {
    for(k=0;k<polygons[p].npt;k++) {
      optimal[p].t[k]=rhs[n];
      double d=fabs(x2[n]-rhs[n]);
      if(d>1.e-06) {
//         printf("node %d: x=%lf x=%lf x=%lf\n",x1[n],x2[n],rhs[n]);
        count++;
        }
      n++;
      }
    }

/*------------------------------------------------------------------------------
  solve for latitudes */  
  n=0;
  for(p=0;p<polygons.size();p++) {
#ifndef PERIODIC
    row=n;
    rhs[row]=y2[n];
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      row=n;
      rhs[row]=K1*y1[n]+K2*y2[n];
      n++;
      }
    row=n;
    rhs[row]=y2[n];
    n++;
#else
    row=n;
    rhs[row]=0;
    n++;
    for(k=1;k<polygons[p].npt;k++) {
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      row=n;
      rhs[row]=K1*y1[n]+K2*y2[n];
      rhs[row]+=K1*y1[n-1]+K1*y1[n+1];
      n++;
      }
#endif
   }
      
  status=LinearSystem_solve(P, rhs);
  
  count=0;
  n=0;
  for(p=0;p<polygons.size();p++) {
    for(k=0;k<polygons[p].npt;k++) {
      optimal[p].p[k]=rhs[n];
      double d=fabs(y2[n]-rhs[n]);
      if(d>1.e-06) {
        count++;
        }
      n++;
      }
    }
    
  status=plg_save("distorded.plg", PLG_FORMAT_SCAN, optimal);

  P.destroy();
  
  delete[] rhs;
  
  delete[] x1;
  delete[] y1;
  
  delete[] x2;
  delete[] y2;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Mesh1DElasticDistorsionT_d(mesh_t & mesh, const vector<plg_t> & polygons, vector<plg_t> & optimal, double* K,  bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int e,k,l,m,n,n1,n2,p,status;
  int row,col;
  size_t pos;
  hypermatrix_t P;
  double *rhs;
  double *x1,*y1,*x2,*y2,*d;
  int count;
  size_t neq, nvtxs, nedges;
  bool force_duplicate=false;
  int verbose=0;
  
  neq=0;
  for(p=0;p<polygons.size();p++) {
    neq+=polygons[p].npt;
    }
  
  for(p=0;p<polygons.size();p++) {
    plg_t tmp;
    optimal.push_back(tmp);
    optimal[p].duplicate(polygons[p]);
    }
  
  nvtxs=neq;
  
  int *vertex=new int[nvtxs];
      
  d=new double[nedges];
  
  double *weight=new double[nedges];
  for(n=0;n<nedges; n++) weight[n]=1.0;

  n=0;
  for(p=0;p<polygons.size();p++) {
    d[n]=0;
    n++;
    for(k=1;k<polygons[p].npt;k++) {
      d[n]=plg_length(polygons[p],k-1,k,0);
      n++;
      }
    }
    
  x1=new double[nvtxs];
  y1=new double[nvtxs];
  
  x2=new double[nvtxs];
  y2=new double[nvtxs];
  
  n=0;
  for(p=0;p<polygons.size();p++) {
    for(k=0;k<polygons[p].npt;k++) {
      x1[n]=polygons[p].t[k];
      y1[n]=polygons[p].p[k];
      n++;
      }
    }
  
  for(n=0;n<nvtxs; n++) vertex[n]=-1;
  
  neq=nvtxs;
  
  P.ordering=new ordering_t;
  
  P.ordering->nrows=neq;
  
  P.ordering->cardinal=new int[neq];
  P.ordering->pointer =new int[neq+1];
  
  n=0;
  for(p=0;p<polygons.size();p++) {
    P.ordering->cardinal[n] =1;
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
      P.ordering->cardinal[n] =3;
      n++;
      }
    P.ordering->cardinal[n] =1;
    n++;
    }
      
  P.ordering->pointer[0]=0;
  for(n=0;n<neq; n++) {
    P.ordering->pointer[n+1]=P.ordering->pointer[n]+P.ordering->cardinal[n];
    }
  
  size_t nnz=P.ordering->pointer[neq];
  P.ordering->incidence=new int[nnz];
  
  n=0;
  for(p=0;p<polygons.size();p++) {
    pos=P.ordering->pointer[n];
    P.ordering->incidence[pos]=n;
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
      pos=P.ordering->pointer[n];
      P.ordering->incidence[pos]  =n-1;
      P.ordering->incidence[pos+1]=n;
      P.ordering->incidence[pos+2]=n+1;
      n++;
      }
    pos=P.ordering->pointer[n];
    P.ordering->incidence[pos]=n;
    n++;
    }
              
  status=matrix_reorder(P.ordering);
   
  P.allocate();
  rhs=new double[neq];
  
  P.assign(0.0);
    
/*------------------------------------------------------------------------------
  set position matrix coefficients */
  n=0;
  for(p=0;p<polygons.size();p++) {
    row=n;
    pos=P.ordering->pos(row,row);
    P.packed[pos] =1.0;
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      row=n;
      e=n-1;
      pos=P.ordering->pos(row,row);
      P.packed[pos] += K[e];
      col=n-1;
      pos=P.ordering->pos(row,col);
      P.packed[pos] -= K[e];
      e=n;
      pos=P.ordering->pos(row,row);
      P.packed[pos] += K[e];
      col=n+1;
      pos=P.ordering->pos(row,col);
      P.packed[pos] -= K[e];
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      pos=P.ordering->pos(row,row);
      P.packed[pos] += weight[n];
      n++;
      }
    row=n;
    pos=P.ordering->pos(row,row);
    P.packed[pos] =1.0;
    n++;
    }
      
  P.distributor=new distribution_t;
  P.distributor->context=SEQUENTIAL_COMPUTING;
      
  status=LinearSystem_initialise(&P, 0, 0, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
      
/*------------------------------------------------------------------------------
  solve for distances */  
  n=0;
  for(p=0;p<polygons.size();p++) {
    row=n;
    rhs[row]=d[n];
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
/*------------------------------------------------------------------------------
      elastic condition with optimal neighbours */
      row=n;
      rhs[row]=weight[n]*d[n];
      n++;
      }
    row=n;
    rhs[row]=d[n];
    n++;
    }
  
  status=LinearSystem_solve(P, rhs);
  
  count=0;
  n=0;
  for(p=0;p<polygons.size();p++) {
    for(k=0;k<polygons[p].npt;k++) {
      optimal[p].t[k]=rhs[n];
      double d=fabs(x2[n]-rhs[n]);
      if(d>1.e-06) {
//         printf("node %d: x=%lf x=%lf x=%lf\n",x1[n],x2[n],rhs[n]);
        count++;
        }
      n++;
      }
    }
    
  status=plg_save("distorded.plg", PLG_FORMAT_SCAN, optimal);

  P.destroy();
  
  delete[] rhs;
  
  delete[] weight;
  
  delete[] d;
  
  delete[] x1;
  delete[] y1;
  
  delete[] x2;
  delete[] y2;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_Mesh1DElasticityT(mesh_t & mesh, double* & K, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
 * 
 * aimed to provide a method for smooth mesh distorsion by computing an elasticity
 * matrix from the original mesh that can be used to recompute vertices positions
 * in the distorded mesh after boundary edges modification
 * 
 * here it is coded for triangle meshes, easily adaptable to any other types of
 * meshes (quadrangles,...)
 * 
 * Next step will be a displacement operator to fit mesh boundaries with a given
 * line or polygon.
 * 
 * Boundary line-wise elasticity to be implemented:
 * 
 *   - 1 fixed vertex per boundary
 *   - equilibrium condition for the others: k1 * L1 + k2 * L2 =0
 * 
 *------------------------------------------------------------------------------ 
  
  Compute boundary edge geometrical elasticity
  --------------------------------------------
  
  Constrained problem:
  
    n boundary vertices equations (to be strictly implemented)

    min potential energy :  E=(sum on neighbours) ki((xn-xi)²+(yn-yi)²]
    
    (x[n],y[n])=r*(x[n-1],y[n-1])+(1-r)*(x[n+1],y[n+1]) ???
    
    It leads to 2 conditions:
           -> dE/dx: (sum on neighbours) ki xn = (sum) ki xi
           -> dE/dy: (sum on neighbours) ki yn = (sum) ki yi
  
    min departure from 1 :  D=(sum on edges) (ki-1)²
  
  
  Unconstrained problem (using Lagrange mutipliers):
    
    min J=(sum on edges) (ki-1)² + (sum on vertices) ln (sum on neighbours j) kj(xn-xj), where ki and ln are unknowns
  
    dJ/dki -> 2(ki-1) + (sum) ln (xn-xi) (2 vertices n involved)
  
    dJ/dli -> (sum on neighbours j) kj(x-xj) = 0 (nngh edges involved)
  
    boundary conditions:
  
    ki = 1 for boundary edges
  
    li = 0 for boundary vertices
  
  
  Reference : https://fr.wikipedia.org/wiki/Multiplicateur_de_Lagrange
  
------------------------------------------------------------------------------*/
{
  int e,k,l,m,n,n1,n2,p,status;
  int row,col;
  size_t pos;
  hypermatrix_t R, P;
  double *rhs, *weight,Energy,*energyP1;
  double *x,*y,*d;
  int count;
  char *flag;
  int nvtxs, nedges;
  size_t neq=0;
  bool force_duplicate=false;
  int verbose=0;
//   bool debug=false;
  
  vector<plg_t> polygons;
  
  flag=new char[mesh.nedges];
  for (n=0;n<mesh.nedges;n++) flag[n]='I';
  
  for (int l=0; l<mesh.nlimits; l++) {
    for (int k=0; k<mesh.limits[l].nedges; k++) {
      n=mesh.limits[l].edges[k];
      flag[n]='M';
      }
    }
    
  status=fe_limits2poly(mesh,polygons,flag,debug);
  
  nvtxs=0;
  for(p=0;p<polygons.size();p++) {
    nvtxs+=polygons[p].npt;
    }
  
  nedges=nvtxs;
  
  int *vertex=new int[nvtxs],*edge=new int[nedges];
  
  double mask=-1;
  
  if(K==0) K=new double[nedges];
    
  x=new double[nvtxs];
  y=new double[nvtxs];
  
  d=new double[nedges];
  
  weight=new double[nedges];
  for(n=0;n<nedges; n++) weight[n]=1.0;

  n=0;
  for(p=0;p<polygons.size();p++) {
    d[n]=0;
    n++;
    for(k=1;k<polygons[p].npt;k++) {
      d[n]=plg_length(polygons[p],k-1,k,0);
      n++;
      }
    }

  n=0;
  for(p=0;p<polygons.size();p++) {
    x[n]=0;
    n++;
    for(k=1;k<polygons[p].npt;k++) {
      x[n]=x[n-1]+d[n];
//       x[n]=polygons[p].t[k];
//       y[n]=polygons[p].p[k];
      n++;
      }
    }
  
  for(n=0;n<nvtxs; n++) vertex[n]=-1;

redo:
/*------------------------------------------------------------------------------
  build rigidity matrix ordering */
  R.ordering=new ordering_t;
  
  neq=nvtxs+nedges;
  R.ordering->nrows=neq;
  
  R.ordering->cardinal=new int[neq];
  R.ordering->pointer =new int[neq+1];
  
  n=0;
  for(p=0;p<polygons.size();p++) {
    R.ordering->cardinal[n]  =1;
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
      R.ordering->cardinal[n]  =2;
      n++;
      }
    R.ordering->cardinal[n]  =1;
    n++;
    }
      
  for(n=0;n<nedges; n++) {
    R.ordering->cardinal[n+nvtxs] = 3;
    }
    
  R.ordering->pointer[0]=0;
  for(n=0;n<neq; n++) {
    R.ordering->pointer[n+1]=R.ordering->pointer[n]+R.ordering->cardinal[n];
    }
  
  size_t nnz=R.ordering->pointer[neq];
  R.ordering->incidence=new int[nnz];
  
/*------------------------------------------------------------------------------
  multipliers lines */
  n=0;
  for(p=0;p<polygons.size();p++) {
/*------------------------------------------------------------------------------
    boundary vertices, self-connection */
    pos=R.ordering->pointer[n];
    R.ordering->incidence[pos]=n;
    n++;
/*------------------------------------------------------------------------------
    interior vertices, edges connected to vertex */
    for(k=1;k<polygons[p].npt-1;k++) {
      pos=R.ordering->pointer[n];
      R.ordering->incidence[pos]=n-1+nvtxs;
      pos=R.ordering->pointer[n];
      R.ordering->incidence[pos+1]=n+nvtxs;
      n++;
      }
/*------------------------------------------------------------------------------
    boundary vertices, self-connection */
    pos=R.ordering->pointer[n];
    R.ordering->incidence[pos]=n;
    n++;
    }
        
/*------------------------------------------------------------------------------
  elastic coefficients lines */
  for(n=0;n<nedges; n++) {
    row=n+nvtxs;
    pos=R.ordering->pointer[row];
    R.ordering->incidence[pos]  = row;
    R.ordering->incidence[pos+1]=n;
    R.ordering->incidence[pos+2]=n+1;
    }
      
  status=matrix_reorder(R.ordering);
  
/*------------------------------------------------------------------------------
  set rigidity matrix and rhs coefficients */
  R.allocate();
  rhs=new double[neq];
  
  R.assign(0.0);
  
/*------------------------------------------------------------------------------
  multipliers lines */
  n=0;
  for(p=0;p<polygons.size();p++) {
    row=n;
    pos=R.ordering->pos(row,row);
    R.packed[pos]=1.0;
    rhs[row]=0;
    n++;
    for(k=1;k<polygons[p].npt-1;k++) {
      double dx;
      row=n;
      m=n-1;
      e=n-1;
      col=e+nvtxs;
      dx=x[n]-x[m];
      pos=R.ordering->pos(row,col);
      R.packed[pos] =dx;
      m=n+1;
      e=n;
      col=e+nvtxs;
      dx=x[n]-x[m];
      pos=R.ordering->pos(row,col);
      R.packed[pos] =dx;
      rhs[row]=0.0;
      n++;
      }
    row=n;
    pos=R.ordering->pos(row,row);
    R.packed[pos]=1.0;
    rhs[row]=0;
    n++;
    }
          
/*------------------------------------------------------------------------------
  elastic coefficients lines */
  for(n=0;n<nedges; n++) {
    row=n+nvtxs;
/*------------------------------------------------------------------------------
    potential energy minimisation */
    n1=n;
    n2=n+1;
    double dx=x[n2]-x[n1];
    col=n1;
    pos=R.ordering->pos(row,col);
    R.packed[pos] = -dx;
    col=n2;
    pos=R.ordering->pos(row,col);
    R.packed[pos] = +dx;
    pos=R.ordering->pos(row,row);
/*------------------------------------------------------------------------------
    departure minimisation */
    R.packed[pos]=weight[n]*2.0;
    rhs[row]=weight[n]*2.0;
    }
    
  R.distributor=new distribution_t;
  R.distributor->context=SEQUENTIAL_COMPUTING;
      
  status=LinearSystem_initialise(&R, 0, 0, SOLVER_ID_UMFPACK, force_duplicate, verbose, debug);
  status=LinearSystem_solve(R, rhs);
  
  count=0;
  for(n=0;n<nedges; n++) {
    row=n+nvtxs;
    K[n]=rhs[row];
    if(K[n]<0.0) {
//       printf("edge %d: K=%lf\n",K[n]);
      weight[n]*=2.0;
      count++;
      }
    else if(K[n]<0.1) {
//       printf("edge %d: K=%lf\n",K[n]);
      weight[n]+=2*(1-K[n]);
      count++;
      }
    }
    
  Energy=0;
  for(n=0;n<nedges; n++) {
    n1=n;
    n2=n+1;
    double dx=x[n2]-x[n1];
//     double dy=y[n2]-y[n1];
//     Energy+=K[n]*(dx*dx+dy*dy);
    Energy+=K[n]*(dx*dx);
    }    

  R.destroy();
  delete[] rhs;
  
  printf("#negative K edges=%d energy=%lf\n",count,Energy);
  
  if(count!=0) goto redo;
  
  delete[] weight;
  
  delete[] d;
  delete[] x;
  delete[] y;
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int eigenvalue_decomposition(double **x, int m, int n, double* & eigenvectors, double* & eigenvalues, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, status;
  
  double mask=+INFINITY;
  
  int ITYPE=1;
  char JOB='V', UPLO='U';
  
  int dim=n;
  
  int LWORK  =1 + 6*dim + 2*dim*dim;
  int LIWORK =3 + 5*dim;
  
  double *WORK,*W;
  int    *IWORK;
  double *A,*B;
  
  WORK =     new double [LWORK];
  IWORK =    new int    [LIWORK];
  W =        new double [dim];
//   A =        new double [dim*dim];
  B =        new double [dim*dim];
  
  for(k=0;k<dim*dim;k++) B[k]=0;
  
  double *C=new double[n*n];
  
  for(l=0;l<n;l++) {
    for(k=0;k<n;k++) {
      covariance_t c=get_covariance(x[l], x[k], mask, m);
      int pos=l*n+k;
      C[pos]=c.covariance;
      }
    }
   
  for(l=0;l<n;l++) {
    int pos=l*n+l;
    B[pos]=1.0;
    }
  
/**----------------------------------------------------------------------------
  then solve the eigenvalue problem */
  dsygvd_(&ITYPE, &JOB, &UPLO, &dim, C, &dim, B, &dim, W, WORK, &LWORK, IWORK, &LIWORK, &status);
  
  delete[] WORK;
  delete[] IWORK;
  delete[] B;
  
  eigenvectors=C;
  eigenvalues=W;

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int pca(mesh_t & mesh, int target, plg_t constraint_raw, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, l, m, status;
  double *eigenvectors[2], *eigenvalues[2];
  double *x,*y,ratio;
  statistic_t s;
  point2D_t p,q;
  double *tmp[2];
  plg_t constraint;
  bool equalize=true;
  char *pj_parameters=0;
  projPJ projection;
  
  double L=plg_minlength(constraint_raw, 0)*1000.0;
//   double L=plg_minlength(constraint_raw, 0)/110;
  
  p=constraint_raw.barycenter();
  projection=assign_StereoOblique(p.x, p.y, pj_parameters);
  
//   printf("\ndefine background working grid: resolution=%lf m, projection=%s\n", L, pj_parameters);

  status=plg_cartesian(projection, &constraint_raw, 1);

  constraint=plg_resample(constraint_raw, L, equalize, debug);
  
  status=plg_spherical(projection, &constraint, 1);
  
  p=constraint.barycenter();
  
  x=new double[constraint.npt];
  y=new double[constraint.npt];
  
  for(k=0;k<constraint.npt; k++) {
    constraint.x[k]=constraint.t[k];
    constraint.y[k]=constraint.p[k];
    x[k]=constraint.x[k]-p.x;
    y[k]=constraint.y[k]-p.y;
    }
    
  for(k=0;k<constraint_raw.npt; k++) {
    constraint_raw.x[k]=constraint_raw.t[k];
    constraint_raw.y[k]=constraint_raw.p[k];
    }
 
  tmp[0]=x;
  tmp[1]=y;
  
/*------------------------------------------------------------------------------
  pca of polygon geometry */
  status=eigenvalue_decomposition(tmp, constraint.npt, 2, eigenvectors[0], eigenvalues[0], debug);
  
  delete[] x;
  delete[] y;
  
  x=new double[mesh.limits[target].nvertex];
  y=new double[mesh.limits[target].nvertex];
  
  for(k=0;k<mesh.limits[target].nvertex; k++) {
    m=mesh.limits[target].vertex[k];
    x[k]=mesh.vertices[m].lon;
    y[k]=mesh.vertices[m].lat;
    }
  
  s=get_statistics(x,mesh.limits[target].nvertex);
  q.x=s.mean;
  
  s=get_statistics(y,mesh.limits[target].nvertex);
  q.y=s.mean;
  
  for(k=0;k<mesh.limits[target].nvertex; k++) {
    m=mesh.limits[target].vertex[k];
    x[k]-=q.x;
    y[k]-=q.y;
    }
    
  tmp[0]=x;
  tmp[1]=y;
  
/*------------------------------------------------------------------------------
  pca of mesh limit geometry */
  status=eigenvalue_decomposition(tmp, mesh.limits[target].nvertex, 2, eigenvectors[1], eigenvalues[1], debug);
  
  for(k=0;k<2;k++) {
    double ratio =sqrt(eigenvalues[0][k]/eigenvalues[1][k]);
    printf("eigenvalues %d: %lf %lf ratio=%lf\n",k,eigenvalues[0][k],eigenvalues[1][k],ratio);
    }
  
/*------------------------------------------------------------------------------
  project mesh limit geometry in its eigen-basis (rotation) */
  for(k=0;k<mesh.limits[target].nvertex; k++) {
    double a,b;
    a=eigenvectors[1][0]*x[k]+eigenvectors[1][1]*y[k];
    b=eigenvectors[1][2]*x[k]+eigenvectors[1][3]*y[k];
    x[k]=a;
    y[k]=b;
    }
    
/*------------------------------------------------------------------------------
  normalize mesh limits to polygon shape factors*/
  for(k=0;k<mesh.limits[target].nvertex; k++) {
    x[k]*=sqrt(eigenvalues[0][0]/eigenvalues[1][0]);
    y[k]*=sqrt(eigenvalues[0][1]/eigenvalues[1][1]);
    }
    
/*------------------------------------------------------------------------------
  back to real world using polygon eigen-basis (rotation) */
  for(k=0;k<mesh.limits[target].nvertex; k++) {
    double a,b;
    a=eigenvectors[0][0]*x[k]+eigenvectors[0][2]*y[k];
    b=eigenvectors[0][1]*x[k]+eigenvectors[0][3]*y[k];
    x[k]=a;
    y[k]=b;
    }
    
  
/*------------------------------------------------------------------------------
  back to real world using polygon barycenter (translation) */
  for(k=0;k<mesh.limits[target].nvertex; k++) {
    m=mesh.limits[target].vertex[k];
    mesh.vertices[m].lon=x[k]+p.x;
    mesh.vertices[m].lat=y[k]+p.y;
    }
  
  
  delete[] x;
  delete[] y;
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int find_attractor(mesh_t & mesh, int target, vector<plg_t> & constraint, vector<int> & exclusion, vector<int> & locked, vector<int> & segment, point2D_t & q, bool & attractor, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

------------------------------------------------------------------------------*/
{
  int count,k,l,m,n,pos,status;
  int m1,m2, s;
  double x,y;
  
  attractor=true;
  
  l=target;
  
  x=0;
  y=0;
  for(k=0;k<mesh.limits[l].nvertex; k++) {
    vector2D_t u;
    m=mesh.limits[l].vertex[k];
    pos=vpos(m,locked);
    if(pos!=-1) {
      attractor=false;
      continue;
//       break;
      }
    s=plg_NearestPolygon(constraint, mesh.vertices[m].lon, mesh.vertices[m].lat, m1, m2, u);
    x+=mesh.vertices[m].lon;
    y+=mesh.vertices[m].lat;
    pos=vpos(s,segment);
    if(pos==-1) segment.push_back(s);
    }
  x/=mesh.limits[l].nvertex;
  y/=mesh.limits[l].nvertex;
  q.x=x;
  q.y=y;

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_distord(mesh_t & mesh, vector<plg_t> & constraint, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  mesh reshaping to align boundary edges positions with constraint polygons
  
  
  https://fr.wikipedia.org/wiki/D%C3%A9composition_en_valeurs_singuli%C3%A8res
  https://fr.wikipedia.org/wiki/Analyse_en_composantes_principales
  https://en.wikipedia.org/wiki/Eigenvalues_and_eigenvectors
  
------------------------------------------------------------------------------*/
{
  int count,k,l,m,n,pos,status;
  int m1,m2, s;
  double d;
  vector2D_t u;
  double *K=0, *Kbdy=0;
  double Kmove, ratio;
  vector<int> locked;
  vector<int> exclusion,exterior,empty;
  point2D_t q;
  bool attractor=true;
  vector<plg_t> polygons, target, optimal;
 
/*------------------------------------------------------------------------------
  compute boundary edge elasticity vector */
  status=fe_Mesh1DElasticityT(mesh, Kbdy, debug);
  
/*------------------------------------------------------------------------------
  compute interior edge elasticity vector */
  status=fe_Mesh2DElasticityT(mesh, K, debug);
  
/*------------------------------------------------------------------------------
  lock position of vertices at split boundaries */
  for(n=0;n<mesh.nedges; n++) {
    if(mesh.edges[n].code==MESH_FLAGGED_EDGE) {
      locked.push_back(mesh.edges[n].extremity[0]);
      locked.push_back(mesh.edges[n].extremity[1]);
      }
    }
  
  l=0;
/*------------------------------------------------------------------------------
  changed 09/03/2016; colateral effects not checked */
//   status=find_attractor(mesh, l, constraint, empty, empty, exclusion, q, attractor, debug);
  status=find_attractor(mesh, l, constraint, empty, locked, exclusion, q, attractor, debug);
  vector<plg_t> tmp;
  for(l=0;l<exclusion.size(); l++) {
    tmp.push_back(constraint[exclusion[l]]);
    }
  status=plg_save("exclusion.plg", PLG_FORMAT_SCAN, tmp);
  
/*------------------------------------------------------------------------------
  island fit specific */
  for(l=1;l<mesh.nlimits; l++) {
    bool attractor=true;
    vector<int> segment;
    double x,y;
    x=0;
    y=0;
    for(k=0;k<mesh.limits[l].nvertex; k++) {
      vector2D_t u;
      m=mesh.limits[l].vertex[k];
      pos=vpos(m,locked);
      if(pos!=-1) {
        attractor=false;
        break;
        }
      s=plg_NearestPolygon(constraint, exclusion, mesh.vertices[m].lon, mesh.vertices[m].lat, m1, m2, u);
      x+=mesh.vertices[m].lon;
      y+=mesh.vertices[m].lat;
      pos=vpos(s,segment);
      if(pos==-1) segment.push_back(s);
      }
    x/=mesh.limits[l].nvertex;
    y/=mesh.limits[l].nvertex;
    q.x=x;
    q.y=y;
    if(segment.size()==1 and attractor) {
/*------------------------------------------------------------------------------
      mesh interior limit fits one polygon */
      point2D_t p=constraint[segment[0]].barycenter();
      double length=1000.*plg_length(constraint[segment[0]], 0);
      printf("limit %d (length=%8.1lf) has polygon attractor %d (length=%8.1lf): lon=%9.4lf lat=%9.4lf\n",
              l, mesh.limits[l].length, segment[0],length, mesh.vertices[m].lon, mesh.vertices[m].lat);
      printf("barycenter lon=%9.4lf lat=%9.4lf lon=%9.4lf lat=%9.4lf\n",
              q.x,q.y,p.x,p.y);
      double ratio=length/mesh.limits[l].length;
//       if(ratio>1.05) continue;
      
      status=pca(mesh, l, constraint[segment[0]], debug);
      
      status= fe_savemesh("mesh-pca.nei",MESH_FILE_FORMAT_TRIGRID, mesh);

//       for(k=0;k<mesh.limits[l].nvertex; k++) {
//         vector2D_t u;
//         m=mesh.limits[l].vertex[k];
//         x=p.x+(mesh.vertices[m].lon-q.x)*ratio;
//         y=p.y+(mesh.vertices[m].lat-q.y)*ratio;
//         mesh.vertices[m].lon=x;
//         mesh.vertices[m].lat=y;
//         locked.push_back(m);
//         }
      }
    }
  
  vector<int> edges;
  for(n=0;n<mesh.nedges; n++) {
    if(mesh.edges[n].code>=1) {
      edges.push_back(n);
      }
    }
  
//   int nsteps=4;
  int nsteps=1;
  for(int step=0; step<nsteps; step++) {

    status=fe_limits2poly(mesh,polygons,0,debug);
 
/*------------------------------------------------------------------------------
  compute nearest constraint polygon intersection, hence new positions */
#pragma omp parallel for private(n,k,s,m,pos)
  for(l=0;l<edges.size(); l++) {
    n=edges[l];
    vector2D_t u;
    for(k=0;k<2;k++) {
      m=mesh.edges[n].extremity[k];
      pos=vpos(m,locked);
      if(pos!=-1) continue;
      s=plg_NearestPolygon(constraint, mesh.vertices[m].lon, mesh.vertices[m].lat, m1, m2, u);
      mesh.vertices[m].lon-=u.x;
      mesh.vertices[m].lat-=u.y;
      }
    }
  
// /*------------------------------------------------------------------------------
//   compute new boundary vertices position */
  status=fe_limits2poly(mesh,target,0,debug);
  status=plg_save("target.plg", PLG_FORMAT_SCAN, target);
// //   Kmove=0.25*(2*step+1.0);
//   Kmove=0.5*(step+1.0);
//   ratio=0.25*(nsteps+3*step+1.0)/(double) nsteps;
//   status=fe_Mesh1DElasticDistorsionT(mesh, polygons, target, optimal, Kbdy, Kmove, ratio, debug);
//   
//   for(l=1;l<mesh.nlimits; l++) {
//     for(k=0;k<mesh.limits[l].nvertex; k++) {
//       n=mesh.limits[l].edges[k];
//       if(mesh.edges[n].code<1) continue;
//       m=mesh.limits[l].vertex[k];
//       pos=vpos(m,locked);
//       if(pos!=-1) continue;
//       mesh.vertices[m].lon=optimal[l].t[k];
//       mesh.vertices[m].lat=optimal[l].p[k];
//       }
//     }
  plg_destroy_entries(polygons);
  polygons.clear();
  plg_destroy_entries(target);
  target.clear();
  plg_destroy_entries(optimal);
  optimal.clear();
  
  }
      
/*------------------------------------------------------------------------------
  compute new interior vertices position */
  status=fe_Mesh2DElasticDistorsionT(mesh, K, debug);
  
  status= fe_savemesh("mesh-distorded.nei",MESH_FILE_FORMAT_TRIGRID, mesh);
  
  delete[] K;
  delete[] Kbdy;
      
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int nndes,status,*elts,fmt,fefmt;
  int i,j,k,l,L,m,n;
  int option,channels=0,smooth=0;
  FILE *file;
  FILE *out;
  char *keyword,*zone=NULL,*gridfile=NULL;
  char *meshfile=NULL,*destruction_polygons=NULL,*output=NULL,*selection_polygons=NULL,*constraint_polygons=NULL;
  char *meshsize=NULL;
  mesh_t mesh,refined,internal,external,*splitted,final;
  double dmax=0,cmax=0,maxrate=0.75;
  int *selected,*targeted,nselected;
  criteria_t criteria;
  plg_t *polygones;
  int npolygones;
  grid_t grid;
  float *density,mask;
  bool debug=false, exact=false;
  int size=1;
  double *depth;
  size_t count;

  vector<plg_t> constraint;

  fct_echo(argc,argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
/* *----------------------------------------------------------------------
    boundary re-sampling*/
    if(strcmp(keyword,"--smooth")==0) {
      smooth=1;
      n++;
      continue;
      }
    if(strcmp(keyword,"--debug")==0) {
      debug=true;
      n++;
      continue;
      }
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          destruction_polygons= strdup(argv[n+1]);
          n++;
          n++;
          break;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          selection_polygons= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          constraint_polygons= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        n++;
        break;
      }
    free(keyword);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  load mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("load mesh ant initialize related tables\n");
  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status!=0) {
      printf("unable to read the original mesh in %s\n",meshfile);
      goto error;
      }
    status=fe_list(&mesh);
    if(status!=0) {
      printf("unable to build the element list from the original mesh\n");
      goto error;
      }
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

  nndes=mesh.nvtxs;
  status= fe_edgetable(&mesh,0,0);
  if(status!=0) {
    printf("unable to build the edge list from the original mesh\n");
    goto error;
    }
    
  status=fe_vertex_crosstables02(&mesh);
  
  status= fe_codetable2(&mesh,0,1,0);
  if(status!=0) {
    printf("unable to rebuild the limits table and codes of the original mesh\n");
    goto error;
    }
    
  if(destruction_polygons!=0) {
    printf("#################################################################\n");
    printf("delete mesh section: %s\n",destruction_polygons);
    status=fe_MeshCut(mesh, destruction_polygons, exact, debug);
    }

  if(exact) {
    status=fe_presplit(mesh, selection_polygons, 0.5, "");
    if(status!=0) goto error;
    size=3;
    }
  else {
    size=1;
    }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  select edges to be refined
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  selected=new int[mesh.nedges];
  targeted=new int[mesh.nedges];
  
  for (n=0;n<mesh.nedges;n++) {
    selected[n]=1;
    }
    
  if(selection_polygons!=0) {
    printf("#################################################################\n");
    printf("select edges from polygons: %s\n",selection_polygons);
    nselected=fe_selectedges_01(mesh, selection_polygons, selected);
    }
  else {
    nselected=mesh.nedges;
    }
  if(nselected==-1) goto error;
  if(nselected== 0) goto end;

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  extract submesh from selection to allow decent cartesian reshape 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("split mesh\n");
  splitted=fe_split(mesh, selected, targeted, 0, size, (string) "mesh-reshape", debug);
  if(splitted==NULL) {
    printf("unable to split the original mesh\n");
    goto error;
    }
    
  printf("#################################################################\n");
  printf("check split connexity\n");
  status=fe_connex(splitted[1]);
  if(status==-1) {
    printf("internal mesh not connex\n");
    goto error;
    }
  
/*------------------------------------------------------------------------------
  transform mesh limits into polygons*/
  status=fe_limits2poly(splitted[1], &polygones, &npolygones,true);
 
  status=plg_load(constraint_polygons, constraint);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  reshape mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_distord(splitted[1], constraint, debug);
      
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  assembly and save final mesh mesh
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  printf("#################################################################\n");
  printf("merge meshes\n");
  final=fe_merge(splitted,(double) 0.01, (char *) 0);
  
  if(debug) 
    status= fe_savemesh("mesh-reshape-no-reshape.nei",MESH_FILE_FORMAT_TRIGRID, final);
  
  if(smooth!=0) {
    printf("#################################################################\n");
    printf("reshape mesh\n");
    status= fe_reshapeall(final,3);
    }

  if(output==0) output=strdup("mesh-reshape.nei");
  status= fe_savemesh(output,MESH_FILE_FORMAT_TRIGRID, final);
  
  depth=new double[final.nvtxs];
  for(n=0;n<final.nvtxs;n++) {
    depth[n]=final.vertices[n].h;
    }
 
  count=occurence(-1.0, depth, mesh.nvtxs);
  
  if(count!=mesh.nvtxs) {
    status=fe_ascii_saver1("depth.s2r", final, depth, -99999.9, LGP1, fmt, (char **) 0);
    }
  
  printf("#################################################################\n");
  printf("summary:\n");
  printf("before upgrade: %d vertices %d elements\n",mesh.nvtxs,mesh.ntriangles);
  printf("after  upgrade: %d vertices %d elements\n",final.nvtxs,final.ntriangles);
  
end: __OUT_BASE_LINE__("end of mesh-reshape ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
