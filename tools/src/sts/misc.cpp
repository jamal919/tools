
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
#include "poc-netcdf.hpp"
#include "netcdf-proto.h"
#include "map.h"
#include "geo.h"
#include <netcdf.h>

extern  int fe_ascii_loadr1(const char *, mesh_t, float *, float *, int);
extern  int fe_ascii_loadr2(const char *, mesh_t, float *, float *, float *, int);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double SquareDst(double lon1, double lat1, double lon2, double lat2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------
  Compute approximate distance
  Longitude and latitude in radians
-----------------------------------------------------------------------------*/
{
  double SquareDstq;
  double dX, dY;

  lon2 = geo_recale(lon2, lon1);//tools update

  dX = R * cos(0.5 * (lat1 + lat2)) * (lon1 - lon2);
  dY = R * (lat1 - lat2);
  SquareDstq = dX * dX + dY * dY;

  return (SquareDstq);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void load_optionals(mesh_t & mesh, state2d_t & state2D, parameter_t & data2D,int fmt, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k, m, n, status;
  int nndes;
  float *h, *N, *r, *dhdx, *dhdy, mask;
  char filename[1024], rootname[1024], msg[1024];

  if(gCPU_ID==gCPU_MASTER) {
    printf(SEPARATOR_1);
    printf("\nLoad optionals input files:\n");
    }

  rootname[0] = 0;

  switch(fmt) {
    case LGP0:
      nndes=mesh.ntriangles;
      break;
    case LGP1:
      nndes=mesh.nvtxs;
      break;
    case NCP1:
      nndes=mesh.nedges;
      break;
    case CQP1:
      nndes=mesh.nvtxs;
      break;
    default:
      check_error(-1, "load_optionals: discretisation not implemented", __LINE__, __FILE__, 1);
      break;
    }


/*----------------------------------------------------------------------
  external file for bottom topography */
  if(strcmp(TopoFile, "NONE") != 0) {
    strcpy(filename, rootname);
    strcat(filename, TopoFile);
    h = new float[nndes];
    if(!h)
      check_error(-1, "allocation failure", __LINE__, __FILE__, 1);
    status = fe_ascii_loadr1(filename, mesh, h, &mask, fmt);
    if(status != 0) {
      sprintf(msg, "topo file input failure : %s", filename);
      check_error(-1, msg, __LINE__, __FILE__, 1);
      }
      
// /*----------------------------------------------------------------------
//     bathymetry offset */
//     status=add_TopographyOffsets(mesh, h);
    
    for(n = 0; n < nndes; n++) {
      data2D.h[n] = h[n];
      }
    zaparr(h);
    
// /*----------------------------------------------------------------------
//     bathymetry scale factor */
//     double scale=tugo_cfg->topography->DepthScaleFactor.numerical_value<double>();
//     for(n = 0; n < gFEmesh[0].nvtxs; n++) {
//       data2D.h[n] *= scale;
//       }
    }
  else {
    for(n = 0; n < gFEmesh[0].nvtxs; n++) {
      data2D.h[n]=gFEmesh[0].vertices[n].h;
      }
    }

/*----------------------------------------------------------------------
  STS: NO bathymetry offset NOR factor */

/*----------------------------------------------------------------------
  bottom rugosity */
  if(strcmp(RugosityFile, "NONE") != 0) {
    strcpy(filename, rootname);
    strcat(filename, RugosityFile);
    r = new float[nndes];
    if(!r)
      check_error(-1, "allocation failure", __LINE__, __FILE__, 1);
    status = fe_ascii_loadr1(filename, mesh, r, &mask, fmt);
    if(status != 0)
      check_error(-1, "rugosity file input failure", __LINE__, __FILE__, 1);
    for(n = 0; n < nndes; n++) {
      data2D.vwr[n] = r[n];
      }
    zaparr(r);
    } 
  else
    for(n = 0; n < nndes; n++) {
      data2D.vwr[n] = 0.;
      }

/*----------------------------------------------------------------------
  average Brunt-Vaisala frequency : given in s^-1, is ok */
  if(strcmp(BVFrequencyFile, "NONE") != 0) {
    strcpy(filename, rootname);
    strcat(filename, BVFrequencyFile);
    N = new float[nndes];
    if(!N)
      check_error(-1, "allocation failure", __LINE__, __FILE__, 1);
    status = fe_ascii_loadr1(filename, mesh, N, &mask, fmt);
    if(status != 0)
      check_error(-1, "file input failure", __LINE__, __FILE__, 1);
    for(n = 0; n < nndes; n++) {
      state2D.N[n] = N[n];
      }
    zaparr(N);
    }
  else {
    for(n = 0; n < nndes; n++) {
//      state2D.N[n] = 2.0e-03;
      state2D.N[n] = BVFrequencyDefault;
      }
    }

/*----------------------------------------------------------------------
  (file) slope is given in m/km, must be divide by 1000. to be really
  dimensionless*/
  if(strcmp(SlopeFile, "NONE") != 0) {
    strcpy(filename, rootname);
    strcat(filename, SlopeFile);
    dhdx = new float[nndes];
    if(!dhdx)
      check_error(-1, "allocation failure", __LINE__, __FILE__, 1);
    dhdy = new float[nndes];
    if(!dhdy)
      check_error(-1, "allocation failure", __LINE__, __FILE__, 1);
    status = fe_ascii_loadr2(filename, mesh, dhdx, dhdy, &mask, fmt);
    if(status != 0)
      check_error(-1, "file input failure", __LINE__, __FILE__, 1);
    for(n = 0; n < nndes; n++) {
//       if(abs(dhdx[n])> 100.) dhdx[n]*=100./abs(dhdx[n]);
//       if(abs(dhdy[n])> 100.) dhdy[n]*=100./abs(dhdy[n]);
      data2D.dhdx[n] = dhdx[n]/1000.;
      data2D.dhdy[n] = dhdy[n]/1000.;
      }
    zaparr(dhdx);
    zaparr(dhdy);
    need_gradient = 0;
    }
  else {
    for(n = 0; n < nndes; n++) {
      data2D.dhdx[n] = 0.;
      data2D.dhdy[n] = 0.;
      }
    need_gradient = 1;
    }

  if(gCPU_ID==gCPU_MASTER) {
    printf("\nLoad optionals input files: OK\n");
    printf(SEPARATOR_2);
    }

}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void node_bw(triangle_t * elt, int nn, int ne, int *hbw, int max_nbr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, ii, iii, j, jjt, jnti, jsub, mem1, swch;
  int *jmem, *memjt;

  jmem  = new int[nn];
  memjt = new int[nn * max_nbr];

/*----------------------------------------------------------------------
    check for consistency
----------------------------------------------------------------------*/

  for(i = 0; i < nn; i++)
    jmem[i] = 0;
  for(i = 0; i < ne; i++) {
    for(j = 0; j < 3; j++) {
      jmem[elt[i].vertex[j]] = 1;
      }
    }
  for(i = 0; i < nn; i++) {
    if(jmem[i] != 1) {
      check_error(-1, "node_bw() found node not in incidence list", __LINE__, __FILE__, 1);
      }
    }

/*----------------------------------------------------------------------
  Establish bandwidth and relationships between nodes
----------------------------------------------------------------------*/

  *hbw = 0;
  for(j = 0; j < nn; j++)
    jmem[j] = 0;

  for(j = 1; j < ne; j++) {
    for(i = 0; i < 3; i++) {
      jnti = elt[j].vertex[i];
      jsub = jnti * max_nbr;
      for(ii = 0; ii < 3; ii++) {
        if(ii == i)
          continue;
        jjt = elt[j].vertex[ii];
        mem1 = jmem[jnti];
        swch = 0;
        for(iii = 0; iii < mem1; iii++) {
          if(memjt[jsub + iii] == jjt) {
            swch = 1;
            break;
          }
        }
        if(swch == 1)
          continue;
        memjt[jsub + jmem[jnti]] = jjt;
        jmem[jnti]++;
        if(abs(jnti - jjt) > *hbw)
          *hbw = abs(jnti - jjt);
      }
    }
  }
  zaparr(jmem);
  zaparr(memjt);
}

// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int element_bw(triangle_t * elt, int nn, int ne, int *hbw, int max_nbr)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int i, ii, n, j, jj, k;
//   int **elist, *nlist, nmax;
//   int **elt_nghl, *elt_nngh, elt_nnghmax;
// 
//   nlist = new int[nn];
// 
// /*----------------------------------------------------------------------
//   Establish element bandwidth 
//   developped to optimize point detection in external mesh
// ----------------------------------------------------------------------*/
// 
//   nmax = 0;
//   for(j = 0; j < nn; j++)
//     nlist[j] = 0;
//   for(j = 0; j < ne; j++) {
//     for(i = 0; i < 3; i++) {
//       n = elt[j].vertex[i];
//       nlist[n]++;
//       nmax = MAX(nmax, nlist[n]);
//       }
//     }
//   printf("maximum elements/node = %d \n", nmax);
// 
//   elist = new int *[nn];
//   for(j = 0; j < nn; j++)
//     elist[j] = new int[nmax];
//   for(j = 0; j < nn; j++)
//     nlist[j] = 0;
// 
//   for(j = 0; j < ne; j++) {
//     for(i = 0; i < 3; i++) {
//       n = elt[j].vertex[i];
//       elist[n][nlist[n]] = j;
//       nlist[n]++;
//       }
//     }
// 
//   *hbw = 0;
// 
//   for(j = 0; j < ne; j++) {
//     for(i = 0; i < 3; i++) {
//       n = elt[j].vertex[i];
//       for(ii = 0; ii < nlist[n]; ii++) {
//         *hbw = MAX(abs(j - elist[n][ii]), *hbw);
//         }
//       }
//     }
// 
//   printf("element half-bandwidth = %d \n", *hbw);
// 
//   elt_nnghmax = 3 * max_nb;
// 
//   elt_nngh = new int[ne];
//   for(j = 0; j < ne; j++)
//     elt_nngh[j] = 0;
//   elt_nghl = new int *[ne];
//   for(j = 0; j < ne; j++)
//     elt_nghl[j] = new int[elt_nnghmax];
// 
//   for(j = 0; j < ne; j++) {
//     for(i = 0; i < 3; i++) {
//       n = elt[j].vertex[i];
//       for(ii = 0; ii < nlist[n]; ii++) {
//         jj = elist[n][ii];
//         for(k = 0; k < elt_nngh[j]; k++)
//           if(elt_nghl[j][k] == jj)
//             goto skip;
//         if(elt_nngh[j] == elt_nnghmax)
//           goto error;
//         elt_nghl[j][elt_nngh[j]] = jj;
//         elt_nngh[j]++;
//       skip:
//         k = 0;
//       }
//     }
//   }
// 
//   elt_nnghmax = 0;
//   for(j = 0; j < ne; j++) {
//     elt_nnghmax = MAX(elt_nngh[j], elt_nnghmax);
//     }
// 
//   printf("element max neighbours = %d (presumedly less than %d) \n",  elt_nnghmax, 3 * max_nb);
// 
//   zaparr(nlist);
//   for(j = 0; j < nn; j++)
//     delete[] elist[j];
//   delete[] elist;
// 
//   for(j = 0; j < ne; j++)
//     delete[] elt_nghl[j];
//   delete[] elt_nghl;
//   delete[] elt_nngh;
// 
//   return (0);
// 
// error:
//   return (-1);
// }

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ShowDimension(FILE * out)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  fprintf(out, "#Number of triangles          : %d \n", gFEmesh[0].ntriangles);
  fprintf(out, "#Number of quadrangles        : %d \n", gFEmesh[0].nquadrangles);
  fprintf(out, "#Number of vertices           : %d \n", gFEmesh[0].nvtxs);
  fprintf(out, "#Number of edges              : %d \n", gFEmesh[0].nedges);

  fprintf(out, "#Elevation boundary nodes     : %d \n", gOpenBCs[0].nnodes);
  fprintf(out, "#Flux boundary nodes          : %d \n", gFluxBCs[0].nnodes);
  
  fprintf(out, "#Number max of neighbours     : %d \n", gFEmesh[0].nnghm);
  fprintf(out, "#Number of sampling points    : %d \n", gSampleSet.n);
  fprintf(out, "#Number of sampling sections  : %d \n", gSectionSet.n);
}

