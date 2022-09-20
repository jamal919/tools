/**************************************************************************
  
  Nonlinear finite element time stepping model 
    
***************************************************************************/

#include "tugo-prototypes.h"

/**@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

  Development notes:

  Check : 

  Note: Mostly obsolete, to be re-worked with hypermatrix ressources


@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@*/
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void gradient_operatorLGP1(const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  int dim;
  int nndes;
  int status;
  float C, hC, hc;
  double S, dP, dQ, dZdx, dZdy;

  nndes = mesh.nvtxs;

/*----------------------------------------------------------------------
  gradient operator initialisation  */
  for(n = 0; n < nndes; n++) {
    for(i = 0; i <= mesh.nnghm; i++) {
      grad_opx[n][i] = 0;
      grad_opy[n][i] = 0;
      }
    }

/*-----------------------------------------------------------------------
  compute contribution of elements weighted by element area */
  for(l = 0; l < mesh.ntriangles; l++) {
    for(j = 0; j < 3; j++) {
      n = mesh.triangles[l].vertex[j];
//      if(n==FE_GHOST_NODE) continue;
      C = gdata2D[0].C[n];
      dP = mesh.triangles[l].DP[j] / 2.0e0; //  no /A, == multiply then divide
      dQ = mesh.triangles[l].DQ[j] / 2.0e0;
      dZdx = +dQ / C;
      dZdy = -dP;
      for(i = 0; i < 3; i++) {
        m = mesh.triangles[l].vertex[i];
//        if(m==FE_GHOST_NODE) continue;
        S = 3. * mesh.vertices[m].mw;   /* total surface weight */
        if(m == n) {
          grad_opx[m][0] += dZdx / S;
          grad_opy[m][0] += dZdy / S;
          }
        else {
          for(k = 0; k < mesh.vertices[m].nngh; k++) {
            if(mesh.vertices[m].ngh[k] == n) {
              grad_opx[m][k + 1] += dZdx / S;
              grad_opy[m][k + 1] += dZdy / S;
              break;
              }
            }
          }
        }
      }
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void gradient_operatorCQP1(const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, j, k, l, m, n;
  int dim;
  int nndes;
  int status;
  float C, hC, hc;
  double S, dP, dQ, dZdx, dZdy;

  nndes = mesh.nvtxs;

/*----------------------------------------------------------------------
  gradient operator initialisation  */
  for(n = 0; n < nndes; n++) {
    for(i = 0; i <= mesh.nnghm; i++) {
      grad_opx[n][i] = 0;
      grad_opy[n][i] = 0;
      }
    }

// /*-----------------------------------------------------------------------
//   compute contribution of elements weighted by element area */
//   for(l = 0; l < mesh.nquadrangles; l++) {
//     for(j = 0; j < 4; j++) {
//       n = mesh.quadrangles[l].vertex[j];
//       C = gdata2D[0].C[n];
//       dP = mesh.triangles[l].DP[j] / 2.0e0; 
//       dQ = mesh.triangles[l].DQ[j] / 2.0e0;
//       dZdx = +dQ / C;
//       dZdy = -dP;
//       for(i = 0; i < 4; i++) {
//         m = mesh.triangles[l].vertex[i];
//         S = 3. * mesh.vertices[m].mw;
//         if(m == n) {
//           grad_opx[m][0] += dZdx / S;
//           grad_opy[m][0] += dZdy / S;
//           }
//         else {
//           for(k = 0; k < mesh.vertices[m].nngh; k++) {
//             if(mesh.vertices[m].ngh[k] == n) {
//               grad_opx[m][k + 1] += dZdx / S;
//               grad_opy[m][k + 1] += dZdy / S;
//               break;
//               }
//             }
//           }
//         }
//       }
//     }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void VertexGradient_operator(const mesh_t & mesh)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    switch(mesh.nature()) {
    case FE_TRIANGLE:
      gradient_operatorLGP1(mesh);
      gP1discretisation=LGP1;
      break;
    case FE_QUADRANGLE:
      gradient_operatorCQP1(mesh);
      gP1discretisation=CQP1;
      break;
    default:
      gradient_operatorLGP1(mesh);
      gP1discretisation=-1;
      break;
    }
}
