
/**************************************************************************

  XSCAN/GENESIS, ocean data editor and viewer, 2006-2009

  Unstructured Ocean Grid initiative (UGO)

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université Dalhousie, Halifax, Canada

E-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

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

#include "legend.def"

#include "tools-structures.h"

#include "legend.h"

#include "geo.h"
#include "map.h"

/*----------------------------------------------------------------------------*/
/* Fonctions appelees depuis Fortran */
/*----------------------------------------------------------------------------*/
#ifdef _add_

extern void lgd_free_ ();
#pragma weak lgd_free_ = lgd_free

#endif

extern int DefaultColumn;

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> legend_t *lgd_import_template(double *x, double *y, T **z, int np, int ncols)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int	k,l,m,n,status,nlgds;
  legend_t   *legends=NULL;
  legend02_t *legend02;
  double d;
  vector<int> start, finish;

  nlgds=1;
  legends=new legend_t[nlgds];

  for(l=0;l<nlgds;l++) {
    legends[l].Type=LGD_GRAPH;
    legends[l].ptr=new char[sizeof(legend02_t)]; /// HERE !!!
    legend02=(legend02_t *) legends[l].ptr;
    legend02->np =np;
    legend02->nz =ncols;
    legend02->points=new point_t[legend02->np];
    m=0;
    for (n=0; n < np;n++) {
      legend02->points[m].N=n;
      legend02->points[m].p=y[n];
      legend02->points[m].t=x[n];
      legend02->points[m].z=new float[legend02->nz];
      for(k=0;k<legend02->nz;k++) legend02->points[m].z[k]=z[k][n];
      m++;
      }
    legend02->ID=l;
    legend02->pen0=0;
    legend02->pen1=2;
    legend02->currentz=0;
    }

  status=LGD_STATUS_OK;
  return (legends);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  legend_t *lgd_import(double *x, double *y, float **z, int np, int ncols)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  legend_t   *legends=NULL;
  legends=lgd_import_template(x, y, z, np, ncols);
  return (legends);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  legend_t *lgd_import(double *x, double *y, double **z, int np, int ncols)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  legend_t   *legends=NULL;
  legends=lgd_import_template(x, y, z, np, ncols);
  return (legends);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  legend_t *lgd_import(const vector<point_t> & points)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int	k,l,m,n,status,nlgds;
  legend_t   *legends=NULL;
  legend02_t *legend02;
  double d;
  vector<int> start, finish;

  nlgds=1;
  legends=new legend_t[nlgds];

  for(l=0;l<nlgds;l++) {
    legends[l].Type=LGD_GRAPH;
    legends[l].ptr=new char[sizeof(legend02_t)]; /// HERE !!!
    legend02=(legend02_t *) legends[l].ptr;
    legend02->np =points.size();
    if(legend02->np>0){
      legend02->nz =points[0].nz;
      legend02->points=new point_t[legend02->np];
      }
    else{
      legend02->nz =-1;
      legend02->points=0;
      }
    m=0;
    for (n=0; n < points.size();n++) {
      legend02->points[m]=points[n];
      m++;
      }
    legend02->ID=l;
    legend02->pen0=0;
    legend02->pen1=2;
    legend02->currentz=0;
    }

  status=LGD_STATUS_OK;
  return (legends);
}


