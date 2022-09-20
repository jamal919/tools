
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2017

  Unstructured Ocean Grid initiative

Developers:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
Contributors:
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France
  Frédéric Dupont    Canada environement, Canada

***************************************************************************/

// #include "tugo-prototypes.h"
// #include "tmatrix.h"
#include "fe.h"
#include "geo.h"
#include <cmath>
#include <list>

#include "constants.h"
#include "functions.h"
// #include "tools-assertions.h"
#include "matrix.h"
#include "periodic.h"
#include "netcdf-proto.h"
// #include "poc-netcdf.h"
#include "poc-netcdf.hpp"


using namespace std;

#define TEST

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int periodic_neighbours(discretisation_t & descriptor, vector<paire_t> paires, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*-----------------------------------------------------------------------------
  neihbour list of periodic paires must be augmented to comply with matrix
  periodization
-----------------------------------------------------------------------------*/
{
  int i, k, n, n1, n2, status;
  int pos;
  
/*-----------------------------------------------------------------------------
  update nodes (periodic nodes neighbour list)*/
  for ( k=0; k< paires.size(); k++ ) {
    n1=paires[k].value[0];
    n2=paires[k].value[1];
    int count2=0;
/**----------------------------------------------------------------------------
    first copy existing n2 neighbour list*/
    int *tab2=new int[descriptor.nodes[n1].nnghbs+descriptor.nodes[n2].nnghbs+2];
    for(i=0;i<descriptor.nodes[n2].nnghbs;i++) {
      tab2[count2]=descriptor.nodes[n2].nghbs[i];
      count2++;
      }
/**----------------------------------------------------------------------------
    then add n1 non-periodized neighbours*/
    for(i=0;i<descriptor.nodes[n1].nnghbs;i++) {
      n=descriptor.nodes[n1].nghbs[i];
/**----------------------------------------------------------------------------
      if mesh is narrow enough, n2 might already be in n1 neighbour list*/
      if(n==n2) continue;
      pos=vpos(n, tab2, count2);
      if(pos==-1) {
        tab2[count2]=n;
        count2++;
        }
      }
    pos=vpos(n1, tab2, count2);
    if(pos==-1) {
      tab2[count2]=n1;
      count2++;
      }
    int count1=0;
/**----------------------------------------------------------------------------
    first copy existing n1 neighbour list*/
    int *tab1=new int[descriptor.nodes[n1].nnghbs+descriptor.nodes[n2].nnghbs+2];
    for(i=0;i<descriptor.nodes[n1].nnghbs;i++) {
      tab1[count1]=descriptor.nodes[n1].nghbs[i];
      count1++;
      }
/**----------------------------------------------------------------------------
    then add n2 non-periodized neighbours*/
    for(i=0;i<descriptor.nodes[n2].nnghbs;i++) {
      n=descriptor.nodes[n2].nghbs[i];
      pos=vpos(n, tab1, count1);
/**----------------------------------------------------------------------------
      if mesh is narrow enough, n2 might already be in n1 neighbour list*/
      if(n==n1) continue;
      if(pos==-1) {
        tab1[count1]=n;
        count1++;
        }
      }
    pos=vpos(n2, tab1, count1);
    if(pos==-1) {
      tab1[count1]=n2;
      count1++;
      }
    if(count1 > 50 || count2 > 50) {
      printf("%s : trouble, n1=%d n2=%d \n",__FUNCTION__,n1,n2);
      }
    delete[] descriptor.nodes[n1].nghbs;
    delete[] descriptor.nodes[n2].nghbs;
    descriptor.nodes[n1].nghbs=tab1;
    descriptor.nodes[n2].nghbs=tab2;
    descriptor.nodes[n1].nnghbs=count1;
    descriptor.nodes[n2].nnghbs=count2;
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int periodic_ordering(ordering_t *ordering, vector<paire_t> paires, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  add periodic connections to matrix ordering

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int i, k, n, n1, n2, status;
  int nrows, pos, nadds=0;
  vptr_t *pointer;
  int *incidence,*cardinal;
  size_t count;
  
  if(paires.size()==0) return(0);
  
  vector<int> adds[2];
  
  nrows     = ordering->nrows;
  
  pointer  = new vptr_t[nrows+1];
  cardinal = new int[nrows];
  
  for(n=0;n<nrows;n++) cardinal[n]=ordering->cardinal[n];

/**----------------------------------------------------------------------------
  scan lines */
  for ( k=0; k< paires.size(); k++ ) {
    n1=paires[k].value[0];
    n2=paires[k].value[1];
    for(i=0; i< ordering->cardinal[n1]; i++) {
      n=ordering->incidence[ordering->pointer[n1]+i];
/**----------------------------------------------------------------------------
      is n already in n2 incidence ? */
      pos=vpos(n, &ordering->incidence[ordering->pointer[n2]], ordering->cardinal[n2]);
      if(pos==-1) {
        adds[1].push_back(n);
        nadds++;
        }
      }
    for(i=0; i< ordering->cardinal[n2]; i++) {
      n=ordering->incidence[ordering->pointer[n2]+i];
/**----------------------------------------------------------------------------
      is n already in n1 incidence ? */
      pos=vpos(n, &ordering->incidence[ordering->pointer[n1]], ordering->cardinal[n1]);
      if(pos==-1) {
        adds[0].push_back(n);
        nadds++;
        }
      }
    cardinal[n1]+=adds[0].size();
    cardinal[n2]+=adds[1].size();
    adds[0].clear();
    adds[1].clear();
    }
  
  if(nadds==0) {
    delete[] pointer;
    delete[] cardinal;
    return(0);
    }
    
  pointer[0]=0;  
  for(n=0;n<nrows;n++) pointer[n+1]=pointer[n]+cardinal[n];
  
  count=ordering->pointer[nrows]+nadds;
  
  incidence=new int[count];
  
  for(n=0;n<count;n++) incidence[n]=-1;
  
  for(n=0;n<nrows;n++) {
    for(i=0; i< ordering->cardinal[n]; i++) {
      incidence[pointer[n]+i]=ordering->incidence[ordering->pointer[n]+i];
      }
    }
  
  for ( k=0; k< paires.size(); k++ ) {
    n1=paires[k].value[0];
    n2=paires[k].value[1];
    for(i=0; i< ordering->cardinal[n1]; i++) {
      n=ordering->incidence[ordering->pointer[n1]+i];
/**----------------------------------------------------------------------------
      is n already in n2 incidence ? */
      pos=vpos(n, &ordering->incidence[ordering->pointer[n2]], ordering->cardinal[n2]);
      if(pos==-1) {
        adds[1].push_back(n);
        }
      }
    for(i=0; i< ordering->cardinal[n2]; i++) {
      n=ordering->incidence[ordering->pointer[n2]+i];
/**----------------------------------------------------------------------------
      is n already in n1 incidence ? */
      pos=vpos(n, &ordering->incidence[ordering->pointer[n1]], ordering->cardinal[n1]);
      if(pos==-1) {
        adds[0].push_back(n);
        }
      }
    for(i=0; i<adds[0].size(); i++) {
      size_t offset=ordering->cardinal[n1];
      incidence[pointer[n1]+i+offset]=adds[0][i];
      }
    for(i=0; i<adds[1].size(); i++) {
      size_t offset=ordering->cardinal[n2];
      incidence[pointer[n2]+i+offset]=adds[1][i];
      }
    adds[0].clear();
    adds[1].clear();
    }

  for(n=0;n<count;n++) {
    if(incidence[n]==-1) {
      printf("%s : anomaly\n",__FUNCTION__);
      }
    }
  for(n=0;n<nrows;n++) {
    for(i=0; i< cardinal[n]; i++) {
      if(incidence[pointer[n]+i]==-1) {
        printf("%s : anomaly %d %d\n",__FUNCTION__,n,i);
        }
      }
    }

  delete[] ordering->pointer;
  delete[] ordering->cardinal;
  delete[] ordering->incidence;
 
  ordering->pointer=pointer;
  ordering->cardinal=cardinal;
  ordering->incidence=incidence;
 
  status=matrix_reorder(ordering);
 
}
  

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int mesh_periodic_geometry(mesh_t *mesh, periodic_t *periodic, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// 
//   extend neighbour list to account for periodic conditions
// 
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   int i, j, k, l, m, mm, m1, m2, n, n1, n2;
//   int alias, count, pos;
//   int edge_id, vertex_id;
//   list<int>::iterator it;
//     
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   
//   create periodic element paire
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   periodic->nelement_paire=periodic->side_1.size();
//   periodic->element_paire=new paire_t[periodic->nelement_paire];
// 
//   for(k=0;k<periodic->side_1.size();k++) {    
//     for(i=0;i<2;i++) {
//       edge_id=periodic->edge_paire[k].value[i];
//       m=mesh->edges[edge_id].shared[0];
//       periodic->element_paire[k].value[i]=m;
//       }
//     }
// 
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   
//   update vertices (periodic vertices neighbour list)
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   for ( k=0; k< periodic->nvertex_paire; k++ ) {
//     n1=periodic->vertex_paire[k].value[0];
//     n2=periodic->vertex_paire[k].value[1];
//     int count2=0;
// /**----------------------------------------------------------------------------
//     first copy existing n2 neighbour list*/
//     int *tab2=new int[mesh->vertices[n1].nngh+mesh->vertices[n2].nngh+2];
//     for(i=0;i<mesh->vertices[n2].nngh;i++) {
//       tab2[count2]=mesh->vertices[n2].ngh[i];
//       count2++;
//       }
// /**----------------------------------------------------------------------------
//     then add n1 non-periodic neighbours*/
//     for(i=0;i<mesh->vertices[n1].nngh;i++) {
//       n=mesh->vertices[n1].ngh[i];
//       if(n==n2) continue;
//       pos=vpos(n,tab2,count2);
//       if(pos==-1) {
//         tab2[count2]=n;
//         count2++;
//         }
//       }
//     pos=vpos(n1,tab2,count2);
//     if(pos==-1) {
//       tab2[count2]=n1;
//       count2++;
//       }
//     int count1=0;
// /**----------------------------------------------------------------------------
//     first copy existing n1 neighbour list*/
//     int *tab1=new int[mesh->vertices[n1].nngh+mesh->vertices[n2].nngh+2];
//     for(i=0;i<mesh->vertices[n1].nngh;i++) {
//       tab1[count1]=mesh->vertices[n1].ngh[i];
//       count1++;
//       }
// /**----------------------------------------------------------------------------
//     then add n2 non-periodic neighbours*/
//     for(i=0;i<mesh->vertices[n2].nngh;i++) {
//       n=mesh->vertices[n2].ngh[i];
//       if(n==n1) continue;
//       pos=vpos(n,tab1,count1);
//       if(pos==-1) {
//         tab1[count1]=n;
//         count1++;
//         }
//       }
//     pos=vpos(n2,tab1,count1);
//     if(pos==-1) {
//       tab1[count1]=n2;
//       count1++;
//       }
//     
//     if(count1>mesh->vertices[n1].nngh+mesh->vertices[n2].nngh+2) TRAP_ERR_EXIT(-1,"anomaly, abort\n");
//     if(count2>mesh->vertices[n1].nngh+mesh->vertices[n2].nngh+2) TRAP_ERR_EXIT(-1,"anomaly, abort\n");
//     
//     delete[] mesh->vertices[n1].ngh;
//     delete[] mesh->vertices[n2].ngh;
//     mesh->vertices[n1].ngh=tab1;
//     mesh->vertices[n2].ngh=tab2;
//     mesh->vertices[n1].nngh=count1;
//     mesh->vertices[n2].nngh=count2;
//     }
//     
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   
//   update edges
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
//   count=0;
//   for ( it=periodic->side_1.begin() ; it != periodic->side_1.end(); it++ ) {
//     k=*it;
//     n1=periodic->edge_paire[count].value[0];
//     n2=periodic->edge_paire[count].value[1];
//     m1=mesh->edges[n1].shared[0];
//     m2=mesh->edges[n2].shared[0];
//     if(m1==m2) {
//       if(debug) printf("%s : periodic edges %d %d have same element %d \n",__FUNCTION__, n1,2,m1);
//       }
//     mesh->edges[n1].shared[1]=m2;
//     mesh->edges[n2].shared[1]=m1;
//     mesh->edges[n1].nshared=2;
//     mesh->edges[n2].nshared=2;
//     mesh->edges[n1].eindex[1]=mesh->edges[n2].eindex[0];
//     mesh->edges[n2].eindex[1]=mesh->edges[n1].eindex[0];
//     count++;
//     }
// 
// /*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
//   
//   force neigbours list's consistency (neighbourship symetry)
//     
//   above vertex manipulation will create non-reciprocal neighbour list, this
//   step is necessary
//   
// xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
// 
// //   debug=true;
// /*------------------------------------------------------------------------------
//   first call to perform the necessary changes */
//   count=fe_chkvertex_consistency(mesh, 1, debug);
//   
// /*------------------------------------------------------------------------------
//   second call to verify all is clear */
//   count=fe_chkvertex_consistency(mesh, 1, debug);
//   
//   if(count!=0) TRAP_ERR_EXIT(-1,"vertex consistency not granted, abort\n");
//   
//   mesh->periodic=1;
//   
//   return(0);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_periodic_LGP0(mesh_t  & mesh, periodic_t *periodic, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, j1, j2, k, l, m, m1, m2, n, n1, n2;
  vector<int> nodes;
  paire_t paire;
    
/*-----------------------------------------------------------------------------
  LGP0 */

  periodic->LGP0_paire.clear();

//   for(k = 0; k < periodic->nedge_paire; k++) {
//     n=periodic->edge_paire[k].value[0];
//     paire.value[0]=mesh.edges[n].shared[0];
//     n=periodic->edge_paire[k].value[1];
//     paire.value[1]=mesh.edges[n].shared[0];
//     periodic->LGP0_paire.push_back(paire);  
//     }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_periodic_LGP1(mesh_t  & mesh, periodic_t *periodic, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, j1, j2, k, l, m, m1, m2, n, n1, n2;
  vector<int> nodes;
  paire_t paire;
    
/*-----------------------------------------------------------------------------
  LGP1 */

  periodic->LGP1_paire.clear();

  for(k = 0; k < periodic->nvertex_paire; k++) {
    periodic->LGP1_paire.push_back(periodic->vertex_paire[k]);
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_periodic_NCP1(mesh_t  & mesh, periodic_t *periodic, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, j1, j2, k, l, m, m1, m2, n, n1, n2;
  vector<int> nodes;
  paire_t paire;
    
/*-----------------------------------------------------------------------------
  NCP1 */

  periodic->NCP1_paire.clear();

  for(k = 0; k < periodic->nedge_paire; k++) {
    periodic->NCP1_paire.push_back(periodic->edge_paire[k]);
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_periodic_CQN1(mesh_t  & mesh, periodic_t *periodic, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, j1, j2, k, l, m, m1, m2, n, n1, n2;
  vector<int> nodes;
  paire_t paire;
    
/*-----------------------------------------------------------------------------
  NCP1 */

  periodic->CQN1_paire.clear();

  for(k = 0; k < periodic->nedge_paire; k++) {
    periodic->CQN1_paire.push_back(periodic->edge_paire[k]);
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mesh_periodic_LGP2(mesh_t  & mesh, periodic_t *periodic, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, j1, j2, k, l, m, m1, m2, n, n1, n2;
  vector<int> nodes;
  paire_t paire;
    
/*-----------------------------------------------------------------------------
  LGP2 */
  periodic->LGP2_paire.clear();
  
  if(mesh.LGP2descriptor.nnodes!=0) {
    for(k=0;k<periodic->nelement_paire;k++) {
      n1=periodic->edge_paire[k].value[0];
      n2=periodic->edge_paire[k].value[1];
      m1=mesh.edges[n1].shared[0];
      m2=mesh.edges[n2].shared[0];
      j1=mesh.edges[n1].eindex[0]+1;
      j2=mesh.edges[n2].eindex[0]+1;
      for(l=0;l<3;l++) {
        i=(2*(j1 % 3) + l) % 6;
        n=mesh.LGP2descriptor.NIbE[m1][i];
        paire.value[0]=n;
        i=(2*(j2 % 3) + 2 - l) % 6;
        n=mesh.LGP2descriptor.NIbE[m2][i];
        paire.value[1]=n;
/*-----------------------------------------------------------------------------
        avoid redundancy due to edge-wise scanning */
        if(vpos(n, nodes) ==-1) {
          nodes.push_back(n);
          periodic->LGP2_paire.push_back(paire);
          }
        }
      }
    }
    
  return(0);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int DiagonalMatrix_periodic_LGP1_template(mesh_t & mesh, periodic_t & periodic, T M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  in case of diagonal matrix, identity condition cannot be implemented

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int alias,j,k,n,row1,row2,col,pos;
  
  if(!periodic.activated) return(0);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  add row2's coefficient to row1's 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  for(k = 0; k < periodic.nvertex_paire; k++) {
    row1=periodic.vertex_paire[k].value[0];
    row2=periodic.vertex_paire[k].value[1];
    M[row1]+=M[row2];
/*------------------------------------------------------------------------------
    for NH solver */
    M[row2]=M[row1]; 
    }
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int DiagonalMatrix_periodic_LGP1(mesh_t & mesh, periodic_t & periodic, double *M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=DiagonalMatrix_periodic_LGP1_template(mesh, periodic, M);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int PackedMatrix_periodic02_template(const mesh_t & mesh, const vector<paire_t> & paires, T & M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  Modify matrix to deal with periodic conditions
  
  Principle: 
  
  in variational approach, matrix coefficents are a summation over elements. For
  a periodised node, the true coefficient is the summation of its own plus the one
  of its alias.
  
  Summation is done for both nodes
  
  WARNING: no RHS modification needed if not the result of element summation
 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int alias,j,k,n,row1,row2,col,pos;
  complex<double> zero=0.0;
  int verbose=0;
  size_t count;
  
  count=0;
  for(k = 0; k < paires.size(); k++) {
    row1=paires[k].value[0];
    row2=paires[k].value[1];
    for(j = M.ordering->pointer[row2]; j < M.ordering->pointer[row2 + 1]; j++) {
      col=M.ordering->incidence[j];
      pos=matrix_absolute_pos(M.ordering, row1, col);
/*------------------------------------------------------------------------------
      periodic mesh possible issue */
      if(M.packed[pos]!=zero && M.packed[j]!=zero) {
        if(verbose==1) printf("%s : trouble at paire %d \n",__FUNCTION__, k);
        count++;
        }
/**----------------------------------------------------------------------------
      add row2's coefficient to row1's */
      M.packed[pos]+=M.packed[j];
/**----------------------------------------------------------------------------
      set row1's coefficient equal to row2's */
      M.packed[j]=M.packed[pos];
      }
    }
  
  if(count!=0) printf("%s : %d paires show anomalies\n",__FUNCTION__, count);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int PackedMatrix_periodic02(const mesh_t &  mesh, const vector<paire_t> & paires, hypermatrix_t & M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
    
  status=PackedMatrix_periodic02_template(mesh, paires, M);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int PackedMatrix_periodic02(const mesh_t &  mesh, const vector<paire_t> & paires, hyperzmatrix_t & M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
    
  status=PackedMatrix_periodic02_template(mesh, paires, M);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int PackedMatrix_periodic01_template(const vector<paire_t> & paires, T & M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  apply periodic condition from 2D paires with 3D implicit expension
   
  Modify matrix to apply periodic conditions from 2D paires with 3D implicit expension
  
  Principle: 
  
  in variational approach, matrix coefficents are a summation over elements. For
  a periodised node, the true coefficient is the summation of its own plus the one
  of its alias.
  
  Summation is done for one node, its alias receive identity condition
  
  WARNING: RHS must be modified accordingly if result of element summation
 
  TODO : implement full multi-variate handling

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int alias,j,k,l,n,p;
  int row2D[2], col2D, row3D[2], col3D, pos;
  size_t hdim,vdim,offset;
  int StopOnError=1;
  
  if(paires.size()==0) return(0);
  
  int rowsize=M.ordering->rowsize;
  int colsize=M.ordering->colsize;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  sum matrix coefficicients in row1 for periodic paires (row1,row2)

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  offset=0;
  for(k=0;k<M.ordering->hdim.size();k++) {
    hdim=M.ordering->hdim[k];
    vdim=M.ordering->vdim[k];
/**----------------------------------------------------------------------------
    add row2's coefficient to row1's */
    for(p = 0; p < paires.size(); p++) {
      row2D[0]=paires[p].value[0];
      row2D[1]=paires[p].value[1];
      for(l=0;l<vdim;l++) {
/*------------------------------------------------------------------------------
        3D expansion */
        row3D[0]=offset+structured3DIndex(row2D[0], l, hdim, vdim);
        row3D[1]=offset+structured3DIndex(row2D[1], l, hdim, vdim);
        pos=0;
        for(j = M.ordering->pointer[row3D[1]]; j < M.ordering->pointer[row3D[1] + 1]; j++) {
          col3D=M.ordering->incidence[j];
/**-----------------------------------------------------------------------
          if M.ordering is monotonically increasing, then acceleration is available*/
//           pos=matrix_absolute_pos(M.ordering, row3D[0], col3D);
//           M.packed[pos]+=M.packed[j];
          pos=matrix_relative_pos2(M.ordering, row3D[0], col3D, pos, 1);
          M.packed[M.ordering->pointer[row3D[0]]+pos]+=M.packed[j];
          }
        int count=0;
        for(j = M.ordering->pointer[row3D[0]]; j < M.ordering->pointer[row3D[0] + 1]; j++) {
          if(M.packed[j]!=0.0) count++;
          }
        if(count==0) {
          printf("matrix has a zeroed line : %d\n",row3D[0]);
          }
        }
      }

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  set identity condition in row2 (0....0 1 0....0 -1 0....0) 

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

    for(p = 0; p < paires.size(); p++) {
      row2D[0]=paires[p].value[0];
      row2D[1]=paires[p].value[1];
      for(l=0;l<vdim;l++) {
/*------------------------------------------------------------------------------
        3D expansion */
        row3D[0]=offset+structured3DIndex(row2D[0], l, hdim, vdim);
        row3D[1]=offset+structured3DIndex(row2D[1], l, hdim, vdim);
        for(j = M.ordering->pointer[row3D[1]]; j < M.ordering->pointer[row3D[1] + 1]; j++) {
          M.packed[j]=0;
          }
        size_t alias=matrix_absolute_pos(M.ordering, row3D[0], row3D[0]);
        pos=matrix_absolute_pos(M.ordering, row3D[1], row3D[1]);
        M.packed[pos]= diagonal;
        pos=matrix_absolute_pos(M.ordering, row3D[1], row3D[0]);
        M.packed[pos]=-diagonal;
        }
      }
    offset+=vdim*hdim;
    }    

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int PackedMatrix_periodic01(const vector<paire_t> & paires, hyperzmatrix_t & M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
    
  status=PackedMatrix_periodic01_template(paires, M);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int PackedMatrix_periodic01(const vector<paire_t> & paires, hypermatrix_t & M)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
    
  status=PackedMatrix_periodic01_template(paires, M);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int PackedMatrix_periodic_template(const periodic_t & periodic, T & M, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  try to implement a minimalist matrix/rhs transformation

  WARNING : in case of multi-variate rhs, works only if all blocks with similar
            discretisation
            
  TODO : implement full multi-variate handling

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status;
  
  vector<paire_t> *paires=0;
  
  if(!periodic.activated) return(0);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  fix mono-variate, mono-dimensional case (hdim and vdim not initialized)
  
  CRITICAL 08/02/2019 : dimension fixing to be handle somewhere else ? 
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  if(M.ordering->vdim.size()==0) {
    M.ordering->vdim.push_back(1);
    }
  if(M.ordering->hdim.size()==0) {
    M.ordering->hdim.push_back(M.ordering->nrows);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  ttach periodic informations; should be held by ordering structure ?
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   M.periodic=&periodic;

  switch(discretisation) {
    case -1:
      printf("%s : matrix %s\n", __FUNCTION__, M.name.c_str());
      check_error(-1, "discretisation field not initialized in hypermatrix_t structure", __LINE__, __FILE__, 1);
      break;
      
    case LGP0:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */ 
      break;
      
    case LGP1:
      M.periodic=&periodic.LGP1_paire;
/*------------------------------------------------------------------------------
      CRITICAL 07/02/2019 : this a pretty bad condition testing, to be quickly revised*/ 
//       if(mesh.nlevels<3) {
//         status=PackedMatrix_periodic_generic(periodic.LGP1_paire, M);
//         }
//       else {
//         status=PackedMatrix_periodic01(mesh, periodic, M);
//         }
      status=PackedMatrix_periodic01(*M.periodic, M);
      break;
      
    case DNP1:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */ 
      break;
      
    case NCP1:
      M.periodic=&periodic.NCP1_paire;
//       status=PackedMatrix_periodic_generic(periodic.NCP1_paire, M);
      status=PackedMatrix_periodic01(*M.periodic, M);
      break;
      
    case LGP2:
      M.periodic=&periodic.LGP2_paire;
//       status=PackedMatrix_periodic_generic(periodic.LGP2_paire, M);
      status=PackedMatrix_periodic01(*M.periodic, M);
      break;
      
    case CQP0:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */ 
//       status=PackedMatrix_periodic_generic(periodic.LGP2_paire, M);
//       status=PackedMatrix_periodic01(periodic.LGP2_paire, M);
      break;
      
    case CQP1:
      TRAP_ERR_EXIT(-1, "%d discretisation not implemented yet", discretisation);
      break;
      
    case CQN1:
      M.periodic=&periodic.CQN1_paire;
//       status=PackedMatrix_periodic_generic(periodic.LGP2_paire, M);
      status=PackedMatrix_periodic01(*M.periodic, M);
      break;
      
    default:
      TRAP_ERR_EXIT(-1, "%d discretisation not implemented yet", discretisation);
      break;
    }
    
#if DEV_SOLVER_ACTIVATE
//   if(paires!=0) {
//     for(int k=0;k<paires->size();k++) {
//       M.periodic.push_back((*paires)[k]);
//       }
//     }
#endif

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int PackedMatrix_periodic(const periodic_t & periodic, hypermatrix_t & M,int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=PackedMatrix_periodic_template(periodic, M, discretisation);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int PackedMatrix_periodic(const periodic_t & periodic, hyperzmatrix_t & M,int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=PackedMatrix_periodic_template(periodic, M, discretisation);

  return(status);
}


// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//  int MomentumSystem_periodic(mesh_t mesh, periodic_t periodic, double **Au, double **Av, double *Bu, double *Bv)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
// 
//     try to implement a minimalist matrix/rhs transformation
// 
// @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
// {
//   int alias,j,k,n,row1,row2,col,pos;
//   double tmp;
//   
//   if(!periodic.activated) return(0);
//   
// /**----------------------------------------------------------------------------
//   add row2's coefficient to row1's */
//   for(k = 0; k < periodic.nvertex_paire; k++) {
//     row1=periodic.vertex_paire[k].value[0];
//     row2=periodic.vertex_paire[k].value[1];
//     for(j=0;j<2;j++) {
//       tmp=Au[row1][j]+Au[row2][j];
//       fwdMatrix_u.x[row1][j]=tmp;
//       fwdMatrix_u.x[row2][j]=tmp;
//       tmp=Av[row1][j]+Av[row2][j];
//       fwdMatrix_v.x[row1][j]=tmp;
//       fwdMatrix_v.x[row2][j]=tmp;
//       }
//     tmp=gRHS_u.x[row1]+gRHS_u.x[row2];
//     gRHS_u.x[row1]=tmp;
//     gRHS_u.x[row2]=tmp;
//     tmp=gRHS_v.x[row1]+gRHS_v.x[row2];
//     gRHS_v.x[row1]=tmp;
//     gRHS_v.x[row2]=tmp;
//     }
//   return(0);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int rhs_periodic2D_template(periodic_t & periodic, T *rhs, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  modify rhs to deal with periodic meshes
  
  try to implement a minimalist matrix/rhs transformation

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,n,row1,row2,status;
  vector<paire_t> *paires;
  
  if(!periodic.activated) return(0);
    
  switch(discretisation) {      
    case LGP0:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */
      return(0);
      break;
      
    case LGP1:
      paires=&(periodic.LGP1_paire);
      break;
      
    case NCP1:
      paires=&(periodic.NCP1_paire);
      break;
      
    case DNP1:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */
      return(0);
      break;
      
    case LGP2:
      paires=&(periodic.LGP2_paire);
      break;
      
    case CQP0:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */
      return(0);
      break;
      
    default:
      TRAP_ERR_EXIT(-1, "%d discretisation not implemented yet", discretisation);
      break;
    }  
  
/**----------------------------------------------------------------------------
  add row2's coefficient to row1's, set identity coefficient in row2 */
  for(k = 0; k < paires->size(); k++) {
    row1=(*paires)[k].value[0];
    row2=(*paires)[k].value[1];
    rhs[row1]+=rhs[row2];
    rhs[row2] =0.0;
    }

// /**----------------------------------------------------------------------------
//   add row2's coefficient to row1's, set identity coefficient in row2 */
//   for(k = 0; k < paires->size(); k++) { /// HERE !!!
//     row1=(*paires)[k].value[0];
//     row2=(*paires)[k].value[1];
//     rhs[row1]+=rhs[row2];
//     rhs[row2] =rhs[row1];
//     }

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int rhs_periodic2D(periodic_t & periodic, double *rhs, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=rhs_periodic2D_template(periodic, rhs, discretisation);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int rhs_periodic2D(periodic_t & periodic, complex<double> *rhs, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=rhs_periodic2D_template(periodic, rhs, discretisation);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int rhs_periodic3D_template(mesh_t & mesh, periodic_t & periodic, T **rhs, int vdim, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  try to implement a minimalist matrix/rhs transformation

  TODO : implement full multi-variate handling

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,n,row1,row2,status;
  vector<paire_t> *paires;
  
  if(!periodic.activated) return(0);
    
  switch(discretisation) {      
    case LGP0:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */
      return(0);
      break;
      
    case LGP1:
      paires=&(periodic.LGP1_paire);
      break;
      
    case NCP1:
      paires=&(periodic.NCP1_paire);
      break;
      
    case DNP1:
/*------------------------------------------------------------------------------
      diagonal (?) */ 
      return(0);
      break;
      
    case LGP2:
      paires=&(periodic.LGP2_paire);
      break;
      
    case CQP0:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */
      return(0);
      break;
      
    default:
      TRAP_ERR_EXIT(-1, "%d discretisation not implemented yet", discretisation);
      break;
    }  
  
/**----------------------------------------------------------------------------
  add row2's coefficient to row1's, set identity coefficient in row2 */
  for(k = 0; k < paires->size(); k++) {
    row1=(*paires)[k].value[0];
    row2=(*paires)[k].value[1];
    for(int l=0; l<vdim; l++) {
      rhs[row1][l]+=rhs[row2][l];
      rhs[row2][l] =rhs[row1][l];
      }
    }

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int rhs_periodic3D(mesh_t & mesh, periodic_t & periodic, complex<double> **rhs, int vdim, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=rhs_periodic3D_template(mesh, periodic, rhs, vdim, discretisation);
  
  return(status);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int rhs_periodic3D(mesh_t & mesh, periodic_t & periodic, complex<double> *rhs, ordering_t *ordering)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int alias,j,k,l,n,row2D[2],row3D[2],col,pos;
//   size_t hdim,vdim;
//   
//   int offset=0;
//   for(j=0;j<ordering->hdim.size();j++) {
//     hdim=ordering->hdim[j];
//     vdim=ordering->vdim[j];
// /**----------------------------------------------------------------------------
//     add row2's coefficient to row1's, set identity coefficient in row2 */
//     for(k = 0; k < periodic.nvertex_paire; k++) {
//       row2D[0]=periodic.vertex_paire[k].value[0]; // only for LGP1 !!!
//       row2D[1]=periodic.vertex_paire[k].value[1];
//       for(l=0;l<vdim;l++) {
//         row3D[0]=offset+structured3DIndex(row2D[0], l, hdim, vdim);
//         row3D[1]=offset+structured3DIndex(row2D[1], l, hdim, vdim);
//         rhs[row3D[0]]+=rhs[row3D[1]];
//         rhs[row3D[1]] =0.0;
//         }
//       }
//     offset+=vdim*hdim;
//     }
//   return(0);
// }

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   template <typename T> int rhs_periodic01_template(periodic_t *periodic, T *rhs, ordering_t *ordering, int discretisation)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// 
//   modify rhs to deal with periodic meshes
//   
//   try to implement a minimalist matrix/rhs transformation
//   
//   WARNING : in case of multi-variate rhs, works only if all blocks with similar
//             discretisation
// 
//   TODO : implement full multi-variate handling
// 
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   int alias,j,k,l,n,row2D[2],row3D[2],col,pos;
//   size_t hdim,vdim;
//   vector<paire_t> *paires;
//   
//   if(periodic==0) return(0);
//   
//   if(!periodic->activated) {
//     TRAP_ERR_EXIT(-1, "periodic pointer set to empty structure");
//     return(0);
//     }
//     
// /*------------------------------------------------------------------------------
//   CRITICAL 08/02/2019 : waiting for a proper explicit periodic handilng*/ 
//   if(discretisation==-1) return(0);
//   
// /*------------------------------------------------------------------------------
//   CRITICAL 08/02/2019 : dimension fixing to be handle somewhere else?*/ 
//   if(ordering->vdim.size()==0) {
//     ordering->vdim.push_back(1);
//     }
//   if(ordering->hdim.size()==0) {
//     ordering->hdim.push_back(ordering->nrows);
//     }
//     
//   
//   switch(discretisation) {      
//     case LGP0:
// /*------------------------------------------------------------------------------
//       discontinuous, nothing to do (?) */
//       return(0);
//       break;
//       
//     case LGP1:
//       paires=&(periodic->LGP1_paire);
//       break;
//       
//     case NCP1:
//       paires=&(periodic->NCP1_paire);
//       break;
//       
//     case DNP1:
// /*------------------------------------------------------------------------------
//       diagonal (?) */ 
//       return(0);
//       break;
//       
//     case LGP2:
//       paires=&(periodic->LGP2_paire);
//       break;
//       
//     case CQP0:
// /*------------------------------------------------------------------------------
//       discontinuous, nothing to do (?) */
//       return(0);
//       break;
//       
//     default:
//       TRAP_ERR_EXIT(-1, "%d discretisation not implemented yet", discretisation);
//       break;
//     }  
//   
//   size_t offset=0;
//   for(j=0;j<ordering->hdim.size();j++) {
//     hdim=ordering->hdim[j];
//     vdim=ordering->vdim[j];
// /**----------------------------------------------------------------------------
//     add row2's coefficient to row1's, set identity coefficient in row2 */
//     for(k = 0; k < (*paires).size(); k++) {
//       row2D[0]=(*paires)[k].value[0];
//       row2D[1]=(*paires)[k].value[1];
//       for(l=0;l<vdim;l++) {
//         row3D[0]=offset+structured3DIndex(row2D[0], l, hdim, vdim);
//         row3D[1]=offset+structured3DIndex(row2D[1], l, hdim, vdim);
//         rhs[row3D[0]]+=rhs[row3D[1]];
//         rhs[row3D[1]] =0.0;
//         }
//       }
//     offset+=vdim*hdim;
//     }
//   return(0);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int rhs_periodic01(periodic_t *periodic, complex<double> *rhs, ordering_t *ordering, int discretisation)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   
//   status=rhs_periodic01_template(periodic, rhs, ordering, discretisation);
//   
//   return(0);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int rhs_periodic01(periodic_t *periodic, double *rhs, ordering_t *ordering, int discretisation)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   
//   status=rhs_periodic01_template(periodic, rhs, ordering, discretisation);
//   
//   return(0);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int rhs_periodic01_template(const vector<paire_t> & paires, T *rhs, ordering_t *ordering, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  modify rhs to deal with periodic meshes
  
  try to implement a minimalist matrix/rhs transformation
  
  WARNING : in case of multi-variate rhs, works only if all blocks with similar
            discretisation

  TODO : implement full multi-variate handling

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int alias,j,k,l,n,row2D[2],row3D[2],col,pos;
  size_t hdim,vdim;
    
/*------------------------------------------------------------------------------
  CRITICAL 08/02/2019 : waiting for a proper explicit periodic handilng*/ 
  if(discretisation==-1) return(0);
  
/*------------------------------------------------------------------------------
  CRITICAL 08/02/2019 : dimension fixing to be handle somewhere else?*/ 
  if(ordering->vdim.size()==0) {
    ordering->vdim.push_back(1);
    }
  if(ordering->hdim.size()==0) {
    ordering->hdim.push_back(ordering->nrows);
    }
    
  
  switch(discretisation) {      
    case LGP0:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */
      return(0);
      break;
      
    case LGP1:
    case NCP1:
      break;
      
    case DNP1:
/*------------------------------------------------------------------------------
      diagonal (?) */ 
      return(0);
      break;
      
    case LGP2:
      break;
      
    case CQP0:
/*------------------------------------------------------------------------------
      discontinuous, nothing to do (?) */
      return(0);
      break;
      
    default:
      TRAP_ERR_EXIT(-1, "%d discretisation not implemented yet", discretisation);
      break;
    }  
  
  size_t offset=0;
  for(j=0;j<ordering->hdim.size();j++) {
    hdim=ordering->hdim[j];
    vdim=ordering->vdim[j];
/**----------------------------------------------------------------------------
    add row2's coefficient to row1's, set identity coefficient in row2 */
    for(k = 0; k < paires.size(); k++) {
      row2D[0]=paires[k].value[0];
      row2D[1]=paires[k].value[1];
      for(l=0;l<vdim;l++) {
        row3D[0]=offset+structured3DIndex(row2D[0], l, hdim, vdim);
        row3D[1]=offset+structured3DIndex(row2D[1], l, hdim, vdim);
        rhs[row3D[0]]+=rhs[row3D[1]];
        rhs[row3D[1]] =0.0;
        }
      }
    offset+=vdim*hdim;
    }
    
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int rhs_periodic01(const vector<paire_t> *paires, complex<double> *rhs, ordering_t *ordering, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  if(paires==0) return(0);
  
  status=rhs_periodic01_template(*paires, rhs, ordering, discretisation);
  
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int rhs_periodic01(const vector<paire_t> *paires, double *rhs, ordering_t *ordering, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  if(paires==0) return(0);
    
  status=rhs_periodic01_template(*paires, rhs, ordering, discretisation);
  
  return(0);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int sum_periodic(mesh_t mesh, periodic_t periodic, double *rhs)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
// 
//   try to implement a minimalist matrix/rhs transformation
//   
//   WARNING : LGP1 only
// 
// cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
// {
//   int alias,j,k,n,row1,row2,col,pos;
//   double tmp;
//   
//   if(!periodic.activated) return(0);
//   
// /**----------------------------------------------------------------------------
//   add row2's coefficient to row1's, set identity coefficient in row2 */
//   for(k = 0; k < periodic.nvertex_paire; k++) {
//     row1=periodic.vertex_paire[k].value[0];
//     row2=periodic.vertex_paire[k].value[1];
//     tmp=rhs[row1]+rhs[row2];
//     rhs[row1]=tmp;
//     rhs[row2]=tmp;
//     }
//   return(0);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename T> int sum_periodic2D_template(mesh_t & mesh, periodic_t & periodic, T *rhs, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  try to implement a minimalist matrix/rhs transformation

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int k,n,row1,row2,status;
  vector<paire_t> *paires;
  
  if(!periodic.activated) return(0);
    
  switch(discretisation) {      
    case LGP0:
/*------------------------------------------------------------------------------
      diagonal (?) */ 
      return(0);
      break;
      
    case LGP1:
      paires=&(periodic.LGP1_paire);
      break;
      
    case NCP1:
      paires=&(periodic.NCP1_paire);
      break;
      
    case DNP1:
      break;
      
    case LGP2:
      paires=&(periodic.LGP2_paire);
      break;
      
    case CQP0:
/*------------------------------------------------------------------------------
      diagonal (?) */ 
      return(0);
      break;
      
    default:
      check_error(-1, "discretisation not implemented yet", __LINE__, __FILE__, 1);
      break;
    }  
  
/**----------------------------------------------------------------------------
  add row2's coefficient to row1's, set identity coefficient in row2 */
  for(k = 0; k < paires->size(); k++) {
    row1=(*paires)[k].value[0];
    row2=(*paires)[k].value[1];
    rhs[row1]+=rhs[row2];
    rhs[row2] =rhs[row1];
    }

  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sum_periodic2D(mesh_t & mesh, periodic_t & periodic, double *rhs, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=sum_periodic2D_template(mesh, periodic, rhs, discretisation);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int sum_periodic2D(mesh_t & mesh, periodic_t & periodic, complex<double> *rhs, int discretisation)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=sum_periodic2D_template(mesh, periodic, rhs, discretisation);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int force_periodic(const mesh_t & mesh, periodic_t & periodic, double *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    try to implement a minimalist matrix/rhs transformation

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int alias,j,k,n,row1,row2,col,pos;
  double tmp;
  
  if(not periodic.activated) return(0);
  
/**----------------------------------------------------------------------------
  add row2's coefficient to row1's, set identity coefficient in row2 */
  for(k = 0; k < periodic.nvertex_paire; k++) {
    row1=periodic.vertex_paire[k].value[0];
    row2=periodic.vertex_paire[k].value[1];
    tmp=rhs[row1];
    rhs[row2]=tmp;
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int check_periodic2D(mesh_t & mesh, vector<paire_t> & paires, complex<double> *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int alias,j,k,n,row1,row2,col,pos;
    
  for(k = 0; k < paires.size(); k++) {
    row1=paires[k].value[0];
    row2=paires[k].value[1];
    complex<double> tmp=rhs[row1]-rhs[row2];
    printf("%d %d %d %lf\n", k, row1, row2, abs(tmp));
    }
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int check_periodic2D(mesh_t & mesh, vector<paire_t> & paires, double *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

    

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int status=0;
  int k, row1, row2;
  size_t count=0;
    
  for(k = 0; k < paires.size(); k++) {
    row1=paires[k].value[0];
    row2=paires[k].value[1];
    double tmp=rhs[row1]-rhs[row2];
    if(tmp!=0) {
      printf("%s : %4d %4d %4d %g %lf %lf\n", __func__, k, row1, row2, abs(tmp), rhs[row1], rhs[row2]);
      count++;
      }
    }
  
  if(count!=0) {
    status=-1;
    }
    
  return(status);
}

