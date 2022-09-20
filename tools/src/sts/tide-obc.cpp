
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

#include  "tugo-prototypes.h"
#include  "poc-netcdf.hpp"
#include "map.h"
#include "tides.h"

#include "legend.def"
#include "legend.h"

extern  int decode_atlasname(char *wave, char **filename);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int TidalBoundaryConditions_02(discretisation_t & descriptor, specification_obsolete_t spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/**----------------------------------------------------------------------------

  Open boundary condition from obs input file

-----------------------------------------------------------------------------*/
{
  int i, j, k, n, n1, n2, exact;
  double d, d1, d2;
  double lat1, lat2, lon1, lon2, lat, lon;
  double amp, phase, tx1, tx2, txi, ty1, ty2, tyi, a1, G1, a2, G2;
  char msg[1024];
  
  if(spec.nnodes==0) return(0);
  
  if(NTidalData==0) {
    sprintf(msg,"tidal OBC expected, but no data available from %s, please check file",TideDataFile);
    check_error(-1, msg, __LINE__, __FILE__, 1);
    }

  for(k=0;k<spec.nnodes;k++) {
    exact = 0;
    d1 = 1.e30;
    d2 = 1.e30;
    n1 = 0;
    n2 = 0;
    n =   spec.nodes[k];
//     lon = gdata2D[0].lon[n];
//     lat = gdata2D[0].lat[n];
/**----------------------------------------------------------------------------
    descriptor position in degrees */
    lon = descriptor.nodes[n].lon*M_PI/180.0;
    lat = descriptor.nodes[n].lat*M_PI/180.0;
    for(i = 0; i < NTidalData; i++) {
      d = SquareDst(lon, lat, TidalData[i].lon,TidalData[i].lat);
      if(d < d1) {
        d2 = d1;
        n2 = n1;
        d1 = d;
        n1 = i;
        }
      else {
        if(d < d2) {
          d2 = d;
          n2 = i;
          }
        }
//      printf("node %d %lf %lf\n",i,sqrt(d1)/1000.0,sqrt(d2)/1000.0);
      }
//    printf("node %d %d %lf %lf\n",k,n,sqrt(d1)/1000.0,sqrt(d2)/1000.0);

    lon1 = TidalData[n1].lon;
    lon2 = TidalData[n2].lon;
    lat1 = TidalData[n1].lat;
    lat2 = TidalData[n2].lat;

    for(j = 0; j < NTidalWave; j++) {
      a1 = (double) TidalData[n1].Ha[j];
      G1 = (double) TidalData[n1].Hg[j];
      a2 = (double) TidalData[n2].Ha[j];
      G2 = (double) TidalData[n2].Hg[j];
      tx1 = a1 * cos(G1);
      ty1 = a1 * sin(G1);
      tx2 = a2 * cos(G2);
      ty2 = a2 * sin(G2);

      txi = (d2 * tx1 + d1 * tx2) / (d1 + d2);
      tyi = (d2 * ty1 + d1 * ty2) / (d1 + d2);

      amp = sqrt(txi * txi + tyi * tyi);
      phase = (amp == 0. ? 0.0 : atan2(tyi, txi));

      if(phase < 0.0)
        phase += pi2;

      spec.tides[k].ztide.a[j] = (float) amp;
      spec.tides[k].ztide.G[j] = (float) phase;

      a1 = (double) TidalData[n1].Ua[j];
      G1 = (double) TidalData[n1].Ug[j];
      a2 = (double) TidalData[n2].Ua[j];
      G2 = (double) TidalData[n2].Ug[j];
      tx1 = a1 * cos(G1);
      ty1 = a1 * sin(G1);
      tx2 = a2 * cos(G2);
      ty2 = a2 * sin(G2);

      txi = (d2 * tx1 + d1 * tx2) / (d1 + d2);
      tyi = (d2 * ty1 + d1 * ty2) / (d1 + d2);

      amp = sqrt(txi * txi + tyi * tyi);
      phase = (amp == 0. ? 0.0 : atan2(tyi, txi));

      if(phase < 0.0)
        phase += pi2;

      spec.tides[k].utide.a[j] = (float) amp;
      spec.tides[k].utide.G[j] = (float) phase;

      a1 = (double) TidalData[n1].Va[j];
      G1 = (double) TidalData[n1].Vg[j];
      a2 = (double) TidalData[n2].Va[j];
      G2 = (double) TidalData[n2].Vg[j];
      tx1 = a1 * cos(G1);
      ty1 = a1 * sin(G1);
      tx2 = a2 * cos(G2);
      ty2 = a2 * sin(G2);

      txi = (d2 * tx1 + d1 * tx2) / (d1 + d2);
      tyi = (d2 * ty1 + d1 * ty2) / (d1 + d2);

      amp = sqrt(txi * txi + tyi * tyi);
      phase = (amp == 0. ? 0.0 : atan2(tyi, txi));

      if(phase < 0.0)
        phase += pi2;

      spec.tides[k].vtide.a[j] = (float) amp;
      spec.tides[k].vtide.G[j] = (float) phase;
      }
    }
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int TidalBoundaryConditions(mesh_t & mesh, discretisation_t & descriptor, specification_obsolete_t & spec)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------

  Open boundary condition: 

----------------------------------------------------------------------*/
{
  int i, j, k, n, n1, n2, status;
  float a, G, a1, G1, a2, G2, x, y;
  tidal_wave wave;
  double zi, zr;
  char filename[1024];
  FILE *out;
  char *sdate;

  if(strcmp(TideDataFile, "NONE") != 0)
    status=TidalBoundaryConditions_02(descriptor, spec);
  else {
    /* STS: NO */
    }

  if(admittance) {
    /* STS: NO */
    }
    
  if(equilibrium) {
    /* STS: NO */
    }

  for(k = 0; k < WaveList.n; k++) {
    sprintf(filename, "%s/%s.processed.obc",gOutputPath, WaveList.waves[k].name);
/**----------------------------------------------------------------------------
    MPI unsafe */
    out = fopen(filename, "w");
    fprintf(out, "%s\n", WaveList.waves[k].name);
    for(i = 0; i < spec.nnodes; i++) {
      n = spec.nodes[i];
      x = descriptor.nodes[n].lon;
      y = descriptor.nodes[n].lat;
      a =spec.tides[i].ztide.a[k];
      G =spec.tides[i].ztide.G[k] * r2d;
      a1=spec.tides[i].utide.a[k];
      G1=spec.tides[i].utide.G[k] * r2d;
      a2=spec.tides[i].vtide.a[k];
      G2=spec.tides[i].vtide.G[k] * r2d;
//       fprintf(out, "%f %f %f %f %f %f %f %f\n", x, y, a, G, a1, G1, a2, G2);
      fprintf(out, "%d %f %f %f %f %f %f %f %f\n", n, x, y, a, G, a1, G1, a2, G2);
      }
    fclose(out);
    }
}
