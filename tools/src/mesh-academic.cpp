

/*******************************************************************************

  T-UGO tools, 2006-2009

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

/*******************************************************************************

  2D mesh generator (academic, simple geometry meshes)

  F. Lyard, 2009, CNRS/LEGOS, Toulouse, France

  e-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <list>

using namespace std;

#include "tools-structures.h"

#include "functions.h"
#include "geo.h"
#include "fe.h"
#include "archive.h"
#include "map.h"
#include "polygones.h"
#include "grd.h"
#include "list.h"
#include "maths.h"
#include "sym-io.h"
#include "map.def"
#include "mass.h"
#include "topo.h"
#include "statistic.h"

#include "zapper.h"     /*  rutin.h contains common utility routines  */

extern  void fe_LGP0_to_LGP1(double *in, double *out, mesh_t mesh);
#include "statistic.h"

#include "academic.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template <typename S, typename T> int topo_synthetic(S synthetic, const char *bathymetry, grid_t grid, grid_t sgrid, T *topo, T topomask, bool flip)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int   i,j,k,l,m,n,status;
  int   ii,jj,imin,imax,jmin,jmax;
  float *buffer,*radius,mask;
  double z,d, maxscale;
  signed char *selected;
  double x,y,dx,dy;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  parse topo instruction
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  synthetic.parse(grid.x[grid.nx-1]-grid.x[0]);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  apply topo instruction
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
/*------------------------------------------------------------------------------
  initialize topo buffer */
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      topo[m]=topomask;
      }
    }

/*------------------------------------------------------------------------------
  interpolate topo */
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      if(not flip) d=grid.x[m]-grid.x[0];
      else d=grid.x[grid.nx-1]-grid.x[m];
      status=synthetic.depth(d,z);
      if(status==0) {
        topo[m]=-z;
        }
      else {
        status=synthetic.depth(d,z);
        }
      }
    }
    
  status=topo_savefield(bathymetry, sgrid, topo, topomask);
 
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
  
  create smoothed topo
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  T *smoothed=new T[sgrid.nx*sgrid.ny];
  float scale=20.0e+03;
  
  sgrid.dx=sgrid.x[1]-sgrid.x[0];
  sgrid.dy=sgrid.y[sgrid.nx]-sgrid.y[0];
  status=map_smooth_latThenLong(sgrid, topo, topomask, scale, smoothed);

  status=topo_savefield("synthetic-smoothed.nc", sgrid, smoothed, topomask);
 
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int resize_grid(const grid_t & cgrid, grid_t & grid, resize_t resize)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int i,j,m,n;
  range_t<double> range;
  
  
  if(resize.neutral) return(0);
    
  grid=cgrid;
  
  grid.nx=resize.end[0]-resize.start[0]+1;
  grid.ny=resize.end[1]-resize.start[1]+1;
  
  grid.x=new double[grid.Hsize()];
  grid.y=new double[grid.Hsize()];
  
  grid.xmin=cgrid.xmin+resize.start[0]*cgrid.dx;
  grid.ymin=cgrid.ymin+resize.start[1]*cgrid.dy;
 
  grid.xmax=grid.xmin+grid.nx*grid.dx;
  grid.ymax=grid.ymin+grid.ny*grid.dy;
 
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      int m=grid.nx*j+i;
      grid.x[m]=grid.xmin+i*grid.dx;
      grid.y[m]=grid.ymin+j*grid.dy;
      }
    }
 
//   status=map_grid_minmax(&grid);
//   sgrid=map_get_spherical(cgrid.projection,cgrid);

  map_printgrid(cgrid);
  map_printgrid(grid);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int trisolverd(int la, double *a, double *b, double *diag, double *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
    /* Local variables */
    int k;
    double tmp, tmpdiag;

/* -------------------------------------------------------------- */
/* TRIAGONAL SYSTEM SOLVER */
/* -------------------------------------------------------------- */
/* first Gauss elimination */

/* At surface */

    /* Function Body */
    tmp = 1. / diag[0];
    a[0] *= tmp;
    rhs[0] *= tmp;

/* at mid depth */

    for (k = 1; k < la-1; k++) {
        tmpdiag = diag[k] - a[k - 1] * b[k];
        tmp = 1. / tmpdiag;
        a[k] *= tmp;
        rhs[k] = (rhs[k] - rhs[k - 1] * b[k]) * tmp;
    }

/* at bottom, solution in rhs */

    tmpdiag = diag[la-1] - a[la-2] * b[la-1];
    tmp = 1. / tmpdiag;
    rhs[la-1] = (rhs[la-1] - rhs[la-2] * b[la-1]) * tmp;
/* -------------------------------------------------------------- */
/* second and final Gauss elimination to surface */
/* solution in rhs */

    for (k = la-2; k >= 0; k--) {
        rhs[k] -= rhs[k + 1] * a[k];
    }
    return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int laplacian_integration_01(int n, double dx, double *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  integrate a numerical second derivative system
  
  2 Dirichlet (value of solution) BCs imposed at both extremities
  
------------------------------------------------------------------------------*/
{
  int status; 
  double *a, *b, *diag;
  double dx2=dx*dx;
  
  a=new double[n];
  b=new double[n];
  
  diag=new double[n];
  for(int i=0;i<n;i++) {
    a[i]=0;
    b[i]=0;
    diag[i]=0;
    }

/*------------------------------------------------------------------------------
  Diriclet line */
  diag[0]  =1.0;
/*------------------------------------------------------------------------------
  Diriclet line */
  diag[n-1]=1.0;

  for(int i=1;i<n-1;i++) {
    a[i]=1.0/dx2;
    b[i]=1.0/dx2;
    diag[i]=-2.0/dx2;
    }
  
  status=trisolverd(n, a, b, diag, rhs);
  
  delete[] a;
  delete[] b;
  delete[] diag;

  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int laplacian_integration_02(int n, double dx, double *rhs)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  integrate a numerical second derivative system
  
  1 Neuman (value of 1st derivative) BCs imposed at first extremity
  1 Dirichlet (value of solution) BCs imposed at last extremity
  
------------------------------------------------------------------------------*/
{
  int status; 
  double *a, *b, *diag;
  double dx2=dx*dx;
  
  a=new double[n];
  b=new double[n];
  
  diag=new double[n];
  for(int i=0;i<n;i++) {
    a[i]=0;
    b[i]=0;
    diag[i]=0;
    }

// /*------------------------------------------------------------------------------
//   Neuman line : a is upper band */
//   diag[0]  =-1.0/dx;
//   a[0]     =+1.0/dx;
//   
// /*------------------------------------------------------------------------------
//   Diriclet line */
//   diag[n-1]=1.0;

/*------------------------------------------------------------------------------
  Neuman line : a is upper band */
  diag[n-1]  =+1.0/dx;
  b[n-1]     =-1.0/dx;
  
/*------------------------------------------------------------------------------
  Diriclet line */
  diag[0]=1.0;

  for(int i=1;i<n-1;i++) {
    a[i]=1.0/dx2;
    b[i]=1.0/dx2;
    diag[i]=-2.0/dx2;
    }
  
  status=trisolverd(n, a, b, diag, rhs);
  
  delete[] a;
  delete[] b;
  delete[] diag;

  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int comodo_depth(grid_t cgrid, float h0, float h1, float *topo, float mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status; 
  int m, n;
  double x, dx,  *rhs;
  FILE *out;
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  COMODO depths prescription given by topo second derivative and boundary values
  
       x < x0  : d2 b/dx2 = 0.
  x0 < x < x1  : d2 b/dx2 = −0.5 ∗ (1 − cos(π(x − x0)/(x1 − x0))) 
  x1 < x < x2  : d2 b/dx2 = −1 + 0.5 ∗ (1 + (x2 − x0)/(x3 − x1)) ∗ (1 + cos(π(x − x2)/(x2 − x1))) 
  x2 < x < x2  : d2 b/dx2 = 0.5(x2 − x0)/(x3 − x1)(1 + cos(π(x − x2)/(x3 − x2))) 
       x > x3  : d2 b/dx2 = 0
  
  x0 = 426km, x1 = 443km, x2 = 479km, x3 = 484km

  depth(x) = 4000. + 3800 ∗ b(x)/ max(b(x))

  suggested method : db/dx is integrated two times (with a ∆x = 25m)
  
  implementation : use implicit solver
  
  dH/dx=(H(i+1)-H(i+1))/dx
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  
  n=cgrid.nx;
  dx=cgrid.dx;
  
  rhs=new double[n];
  double x0,x1,x2,x3;
  
  x0=cgrid.x[0]+426000;
  x1=cgrid.x[0]+443000;
  x2=cgrid.x[0]+479000;
  x3=cgrid.x[0]+484000;

  //printf("%lf %lf\n", cgrid.x[0], cgrid.x[1]);
  for(int i=0; i<n; i++) {
    m=i;
    x=cgrid.x[m];
    rhs[i]=0.0;
    if (x < x0){
      rhs[i]=0;
      }
    else if (x < x1){
      rhs[i]=-0.5*(1-cos(M_PI*(x-x0)/(x1-x0)));
      }
    else if (x < x2){
      rhs[i]=-1+0.5*(1+(x2-x0)/(x3-x1))*(1+cos(M_PI*(x-x2)/(x2-x1)));
      }
    else if(x < x3){
      rhs[i]=0.5*(x2-x0)/(x3-x1)*(1+cos(M_PI*(x-x2)/(x3-x2)));
      }
    else {
      rhs[i]=0;
      }
    }
    
  for(int i=0; i<n; i++)  rhs[i]*=-1.0;
  
//   rhs[0]  =h0;
//   rhs[n-1]=h1;
//   
//   status=laplacian_integration_01(n, dx, rhs);
  
  rhs[0]  =0;
  rhs[n-1]=0;
  
  status=laplacian_integration_02(n, dx, rhs);
  
  range_t<double> r=poc_minmax(rhs, n, NAN);
  
  for(int i=0; i<n; i++) {
    rhs[i]=4000.-3800.*rhs[i]/r.max;
    }
    
  for(int j=0; j<cgrid.ny; j++) {
    for(int i=0; i<cgrid.nx; i++) {
      m=cgrid.nx*j+i;
      topo[m]=-rhs[i];
      }
    }
    
  out=fopen("comodo-depth.gnu","w");
  for(int i=0; i<n; i++) {
    fprintf(out, "%lf %lf\n", cgrid.x[i], -rhs[i]);
    }
  fclose(out);
  
  delete[] rhs;

  return 0;
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  plg_t plg_sample(plg_t target, grid_t & grid, float *density, float mask, criteria_t & criteria, double tolerance, int equalize, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_t *tmp,out;
  int concat[2];
  double km2m=1000.0;
  double radius=criteria.cellsize*km2m;
  
  equalize=1;
  
  tmp=new plg_t[target.npt-1];
  
  concat[0]=0;
  
  for(size_t k=0;k<target.npt-1;k++) {
    
    if(target.flag[k]=='M') {
      
/*------------------------------------------------------------------------------
      next step necessary to insure new boundary to lon/lat aligned */
      double check=target.cartesian_size(k);
      int n=check/1000.+2;
      plg_t p=plg_t(target.t[k],target.p[k],target.t[k+1],target.p[k+1],n,PLG_SPHERICAL);
      
      status=plg_cartesian(grid.projection, &p, 1);
      status=plg_spherical(grid.projection, &p, 1);
    
// /*------------------------------------------------------------------------------
//       compute min radius */
//       radius=1.e+10;
//       for(size_t l=0;l<p.npt-1;l++) {
//         double C=plg_resolution(p, l, grid, density, mask, criteria,  equalize);
//         radius=min(C,radius);
//         }

      double C=plg_resolution(p, 0, grid, density, mask, criteria,  equalize);
      vector<double> L;
      double L1;
      plg_t q;
      
      equalize=0;
      
      int count=0;
      while(true) {
        L.push_back(C);
        q=plg_resample(p, L, equalize, debug);
        C=plg_resolution(q, count, grid, density, mask, criteria,  equalize);
        while(C < L[count]) {
          L[count]*=0.99;
          q.destroy();
          q=plg_resample(p, L, equalize, debug);
          C=plg_resolution(q, count, grid, density, mask, criteria,  equalize);
          }
        count++;
        if(count==q.npt-1) break;
        C=plg_resolution(q, count, grid, density, mask, criteria,  equalize);
        q.destroy();
        }
      
      equalize=1;
      q.destroy();
      q=plg_resample(p, L, equalize, debug);
      status=plg_spherical(grid.projection, &q, 1);
    
/*------------------------------------------------------------------------------
      equalize sampling */
      for(int pass=0;pass<10;pass++) 
        for(size_t l=1;l<q.npt-1;l++) {
          double L1=q.cartesian_size(l-1);
          double L2=q.cartesian_size(l);
          double L=MAX(L1,L2);
          double C1=plg_resolution(q, l-1, grid, density, mask, criteria,  equalize);
          double C2=plg_resolution(q, l,   grid, density, mask, criteria,  equalize);
          double C=C1+C2, r=1.-C1/C;
          q.x[l]=r*q.x[l-1]+(1.-r)*q.x[l+1];
          q.y[l]=r*q.y[l-1]+(1.-r)*q.y[l+1];
          }
      status=plg_spherical(grid.projection, &q, 1);
      q.SetFlag('M');
      
      tmp[k].duplicate(q);
      if(k!=0) {
        concat[1]=k;
/*------------------------------------------------------------------------------
        aggregate new pieces in first polygon*/
        status=plg_concat(concat, tmp, k+1, PLG_CARTESIAN);
        if(status!=0) {
          printf("concatenation issue : %d\n",k);
          }
        }
      p.destroy();
      q.destroy();

      }
    else {
      plg_t p=plg_t(target.x[k],target.y[k],target.x[k+1],target.y[k+1]);
      p.SetFlag('T');
      tmp[k].duplicate(p);
      if(k!=0) {
        concat[1]=k;
        status=plg_concat (concat, tmp, k+1, PLG_CARTESIAN);
        if(status!=0) {
          printf("concatenation issue : %d\n",k);
          }
        }
      p.destroy();      
      }
    }
  
  out.duplicate(tmp[0]);
  
  plg_deletep(&tmp,target.npt-1);
  
  return(out);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_create_boundaries(plg_array_t *boundaries, plg_array_t plg, double radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int k;
  int equalize=0;

  boundaries->p=new plg_t[plg.n];
  for(k=0;k<plg.n;k++) {
    boundaries->p[k]=plg_resample(plg.p[k],radius, equalize);
    }
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int fe_create_interiors(plg_array_t boundaries, double radius)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  __CHKERR_LINE_FILE__(ENOEXEC,"not finished");exit(ENOEXEC);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_create_rectangle(plg_array_t *plg, double t0, double p0, double dt, double dp)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  int i,j,k;
  double x[4],y[4];
//   double L,H,dt,dp,t0,p0,radius;
  static double cycle[4][2]={(0,0),(1,0),(1,1),(0,1)};

  for(k=0;k<4;k++) {
    x[k]=t0+cycle[k][0]*dt;
    y[k]=p0+cycle[k][1]*dp;
    }

  plg->n=4;
  plg->p=new plg_t[plg->n];

  for(k=0;k<4;k++) {
    i=k;
    j=(k+1)%4;
    plg->p[k]=plg_t(x[i],y[i],x[j],y[j]);
    }

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int plg_create_circle(plg_array_t *plg, double radius, double t0, double p0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int i,j,k,utm_zone,ntoken;
  double *x,*y,x0,y0;
  double L,H;
  char *proj4_options;
  const char *hemisphere;
  projPJ proj;
  char **parms;
  
  L=2.*M_PI*radius;

  plg->n=1;
  plg->p=new plg_t[plg->n];

  plg->p[0].npt=100;
  
  x=new double[plg->p[0].npt];
  y=new double[plg->p[0].npt];

  plg->p[0].t=new double[plg->p[0].npt];
  plg->p[0].p=new double[plg->p[0].npt];
  
  if(p0<0) hemisphere="+south";
  else hemisphere="";
  
  utm_zone=(t0+180)/6+1;
  
  asprintf(&proj4_options,"+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs %s",utm_zone,hemisphere);

  proj = pj_init_plus(proj4_options);
  if (!proj) __TRAP_ERR_EXIT__(1,"Projection initialization failed\n");
  
  geo_to_projection(proj, p0, t0, &x0, &y0);
  pj_free(proj);
  
  for(n=0;n<plg->p[0].npt;n++) {
    double alpha=n*2*M_PI/(plg->p[0].npt-1.);
    x[n]=x0+radius*cos(alpha);
    y[n]=y0+radius*sin(alpha);
    }
    
  status=projection_to_geo (proj4_options, x, y, plg->p[0].npt);
  
  free(proj4_options);
  
  for(n=0;n<plg->p[0].npt;n++) {
    plg->p[0].t[n]=x[n];
    plg->p[0].p[n]=y[n];
    }

  status=plg_save("circle.plg",PLG_FORMAT_SCAN,plg->p,1);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_create_polar(plg_t & p, double radius, double t0, double p0, int nsectors)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int i,j,k,zone,ntoken;
  double *x,*y,x0,y0;
  double L,H;
  char *proj4_options;
  const char *hemisphere;
  projPJ proj;
  char *pj_parameters=new char[256];
  
  L=2.*M_PI*radius;

  p.init(nsectors+1, PLG_INIT_SEPARATE);
  
  x=new double[p.npt];
  y=new double[p.npt];
  
  if(p0<0) hemisphere="+south";
  else hemisphere="";
  
  zone=(t0+180)/3;
  
  proj = assign_StereoNorth(p0, 90.0, t0, pj_parameters);
  
  if (!proj) __TRAP_ERR_EXIT__(1,"Projection initialization failed\n");
  
  geo_to_projection(proj, p0, t0, &x0, &y0);
  
  for(n=0;n<p.npt;n++) {
    double alpha=n*2*M_PI/(p.npt-1.);
    x[n]=x0+radius*cos(alpha);
    y[n]=x0+radius*sin(alpha);
    }
    
  for(n=0;n<p.npt;n++) {
    p.x[n]=x[n];
    p.y[n]=y[n];
    }
    
  status=projection_to_geo (proj, x, y, p.npt);
  
//   free(proj);
  pj_free(proj);

  for(n=0;n<p.npt;n++) {
    p.t[n]=x[n];
    p.p[n]=y[n];
    }
    
//   status=plg_save("circle.plg",PLG_FORMAT_SCAN, &p,1);

  return(status);}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_create_shape(plg_array_t *plg)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  double dt=0.0,dp=0.0,t0=0.0,p0=0.0;
  
  status=plg_create_rectangle(plg, t0,p0,dt,dp);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_create_circularT(bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  aimed to study north pole issues

------------------------------------------------------------------------------*/
{
  int i,j,k,l,m,n,status;
  vector<plg_t> plg, limits;
  plg_t *p;
  int nsectors=16;
  double lon, radius;
  
  char *proj4_options;
  projPJ proj;
  char *pj_parameters=new char[256];
  
  mesh_t *mesh, work;
  
  p=new plg_t;
  nsectors=24;
  lon=180.0/nsectors;
  radius=12500;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
  p=new plg_t;
//   nsectors=16;
  lon=0;
  radius+=7000;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
  p=new plg_t;
//   nsectors=16;
  lon=180.0/nsectors;
  radius+=10000;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
  p=new plg_t;
//   nsectors=16;
  lon=0;
  radius+=12500;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
//   goto finish;
  
  p=new plg_t;
//   nsectors=16;
  lon=180.0/nsectors;
  radius+=15000;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
  p=new plg_t;
  nsectors*=2;
  lon=0.0;
  radius+=15000;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
  p=new plg_t;
//   nsectors=32;
  lon=180.0/nsectors;
  radius+=20000;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
  p=new plg_t;
//   nsectors=32;
  lon=0.0;
  radius+=25000;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
  p=new plg_t;
//   nsectors=32;
  lon=180.0/nsectors;
  radius+=25000;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
  p=new plg_t;
//   nsectors=32;
  lon=0.0;
  radius+=25000;
  status= plg_create_polar (*p, radius, lon, +90.0, nsectors);
  plg.push_back(*p);
  
finish:
  status=plg_save("circle.plg",PLG_FORMAT_SCAN, plg);
  
  
  proj = assign_StereoNorth(90.0, 90.0, 0.0, pj_parameters);
  if (!proj) __TRAP_ERR_EXIT__(1,"Projection initialization failed\n");
  
  status=plg_cartesian(proj,plg);
  
  size_t count=1;
  for(int s=0;s<plg.size()-1;s++) {
    count+=plg[s].npt-1;
    }
    
  point_t* points=new point_t[count];

  count=0;
  points[count].t=0.0;
  points[count].p=0.0;
  count++;
  
  for(int s=0;s<plg.size()-1;s++) {
    for(int k=0;k<plg[s].npt-1;k++) {
      points[count].t=plg[s].x[k];
      points[count].p=plg[s].y[k];
      count++;
      }
    }
    
//   plg[plg.size()-1].npt--;
  limits.push_back(plg[plg.size()-1]);

  status=fe_createnodes(limits, points, count, work);
  delete[] points;
  
  mesh=fe_triangulate(work, proj, 0, debug);
  work.destroy();
  
  for (size_t n=0;n<mesh->nvtxs;n++) {
    double x=mesh->vertices[n].lon;
    double y=mesh->vertices[n].lat;
    projection_to_geo(proj, &(mesh->vertices[n].lat), &(mesh->vertices[n].lon), x, y);
    } 
  status=fe_geometry(mesh);

  status=fe_edgetable(mesh,0,0);
  status=fe_vertex_crosstables02(mesh);
  status=fe_codetable2(mesh,0,1,1);
  
  status=fe_savemesh("mesh-academic.nei",MESH_FILE_FORMAT_TRIGRID, *mesh);
  
  int n1=0;
  int n2=0;
  status=fe_setOBCflags(*mesh, n1, n2, MESH_ELEVATION_NODE);
  status=fe_write_boundarycode("mesh-academic.bel", *mesh, 1);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_create_rectangularT(criteria_t criteria, mesh_t & mesh, plg_t *limits, int nlimits, char *meshsize, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,k,l,m,n,status;
  int count,s;
  mesh_t *spherical;
  double L,H,dt,dp,t0,p0,radius;
  plg_array_t plg,boundaries;
  point_t *points;
  float *density,mask;
  grid_t grid;
  plg_t *polygones=0;
  int npolygones=0;
  const char *output="density.nc";
//   bool debug=true;
  
// /* *----------------------------------------------------------------------------
//   create a geometry (sort of synthetic "shorelines") */
//   status=plg_create_shape(&plg, shape);
  
/**-----------------------------------------------------------------------------
  echo criteria file */
  status= fe_save_criteria("echo.crt", criteria);
  status= fe_save_criteria("mesh-academic-echo.crt", criteria);

/* *-----------------------------------------------------------------------------
  mesh size mapping*/
  printf("#################################################################\n");
  printf("mesh resolution setup\n");
  if(meshsize==0) {
    printf("compute mesh density\n");
    status=fe_ComputeMeshsize(output, criteria, &grid, &density, &mask, limits, nlimits,1, false);
    }
  else {
    printf("load mesh density from %s\n",meshsize);
    status=fe_ReloadMeshsize(meshsize, &grid, &density, &mask, limits, nlimits);
    status= defproj(&grid, limits, nlimits);
    }
    
  printf("#################################################################\n");
  printf("mesh limits sampling\n");
  
  for(int s=0; s<  nlimits; s++) {
    if(limits[s].npt==0) {
      printf("empty polygon %d\n",s);
      continue;
      }
    limits[s].flag=new char[limits[s].npt-1];
    for(int k=0;k<limits[s].npt-1;k++) limits[s].flag[k]='M';

    plg_t *p=new plg_t[limits[s].npt-1];
    plg_t *q=new plg_t[limits[s].npt-1];
    for(int k=0;k<limits[s].npt-1;k++) {
      p[k].duplicate(limits[s],k,k+1);
/* *-----------------------------------------------------------------------------
      necessary to insure mesh periodicity, to be improved*/
      if(k==2) status=plg_flip(p[k]);
      q[k]=plg_sample(p[k], grid, density, mask, criteria, 0.95, 1, debug);
      if(q[k].npt<2) {
        printf("%s : sampling failed\n",__func__);
        printf("deep  size=%lf %lf km\n", criteria.maxsize, criteria.minsize);
        printf("shelf size=%lf %lf km \n", criteria.shelf_maxsize, criteria.shelf_minsize);
        status=plg_print(limits[s]);
        status=plg_print(p[k]);
        return(-1);
        }
      }
    
/* *-----------------------------------------------------------------------------
    alternatively this should be done by applying a translation*/
//     int N=q[2].npt-1;
//     for(size_t l=1;l<N;l++) {
//       q[2].t[l]=q[0].t[N-l]; //lon
//       q[2].x[l]=q[0].x[N-l]; // x
//       }

    for(int k=1;k<limits[s].npt-1;k++) plg_concat(q[0],q[k]);

    limits[s].duplicate(q[0]);
    }
  status=plg_save("mesh-academic.adjusted.plg", PLG_FORMAT_SCAN,  limits, nlimits);

// /* *----------------------------------------------------------------------------
//   create mesh boundaries on the geometry */
//   status=fe_create_boundaries(&boundaries, plg, radius);

/* *-----------------------------------------------------------------------------
  mesh generation*/
  printf("#################################################################\n");
  printf("mesh generation\n");
  mesh=fe_nodit(criteria, grid, density, mask, limits, nlimits, debug);
  status= fe_savemesh("mesh-academic-no-reshape.nei",MESH_FILE_FORMAT_TRIGRID, mesh);

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename TYPE> int fe_create_rectangularQ(grid_t grid,grid_t topogrid, TYPE *topo, TYPE topomask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  size_t count=0;
  char *landmask;
  mesh_t mesh;
  double t,p;
  
  landmask=new char[grid.nx*grid.ny];
 
  for (size_t j=0;j<grid.ny;j++) {
    for (size_t i=0;i<grid.nx; i++) {
      float h;
      double x=map_grid_x(grid,i,j);
      double y=map_grid_y(grid,i,j);
      x=map_recale(topogrid,x);
      status=map_interpolation(topogrid,topo,topomask,x,y,&h);
      size_t m=grid.nx*j+i;
      if(h>0) {
        landmask[m]=0;
        }
      else {
        landmask[m]=1;
        count++;
        }
      }
    }

  status=fe_FGrid2Quadrangle(grid, landmask, topo,topomask,mesh,"quadrangle.nc",debug);

  status=fe_vertex_crosstables02(&mesh);

  int stopon_EdgeError =1;
  int stopon_PinchError=1;
  status=fe_codetable1(&mesh, 0, stopon_EdgeError, stopon_PinchError);

  status=fe_savemesh("test.nei",  MESH_FILE_FORMAT_TRIGRID, mesh);
  
  t=+5.;
  p=-2.;
  int n1=fe_nearest_vertex (mesh,t,p,0,0);

  t=+5.;
  p=+2.;
  int n2=fe_nearest_vertex (mesh,t,p,0,0);
  
  status=fe_setOBCflags(mesh, n1, n2, MESH_ELEVATION_NODE);
  status=fe_write_boundarycode("test.bel", mesh, 1);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_create_hexagonal(grid_t sgrid, grid_t cgrid, geo_t projection, criteria_t criteria, double dx, double length, int y_size, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------
  
  create rectangular shape, fix resolution, hexagonal distribution mesh
    
------------------------------------------------------------------------------*/  
{
  int i,j,k,m,n,status;
  mesh_t mesh;
  double x0,y0,x,y,d;
  
//   double dx=criteria.shelf_minsize;
  
  i=1;
  j=0;
  m=cgrid.nx*j+i;
  x0=cgrid.x[m];
  
  i=1;
  j=cgrid.ny/2;
  m=cgrid.nx*j+i;
  y0=cgrid.y[m];
//   
//   i=cgrid.nx-2;;
//   j=0;
//   m=cgrid.nx*j+i;
//   length=cgrid.x[m]-x0;
//   
  int x_size=1+length/dx;
  
  dx=length/(x_size-1);
  
  
  double dy=dx*sqrt(3.0)/2.0;
  
  int half=y_size/2;
  n=x_size+2*half*x_size;
  
  half=(y_size+1)/2;
  n+=+2*half*(x_size+1);
  
  mesh.vertices=new vertex_t[n];
  
  n=0;
  for(i=0;i<x_size;i++) {
    mesh.vertices[n].lon=x0+i*dx;
    mesh.vertices[n].lat=y0;
    n++;
    }
    
  half=y_size/2;  
  for(j=0;j<half;j++) {
    y=y0-2*(j+1)*dy;
    for(i=0;i<x_size;i++) {
      mesh.vertices[n].lon=x0+i*dx;
      mesh.vertices[n].lat=y;
      n++;
      }
    y=y0+2*(j+1)*dy;
    for(i=0;i<x_size;i++) {
      mesh.vertices[n].lon=x0+i*dx;
      mesh.vertices[n].lat=y;
      n++;
      }
    }
  
  half=(y_size+1)/2;  
  for(j=0;j<half;j++) {
    y=y0-(2*j+1)*dy;
    for(i=0;i<x_size-1;i++) {
      mesh.vertices[n].lon=x0+i*dx+dx/2.0;
      mesh.vertices[n].lat=y;
      n++;
      }
    mesh.vertices[n].lon=x0;
    mesh.vertices[n].lat=y;
    n++;
    mesh.vertices[n].lon=x0+length;
    mesh.vertices[n].lat=y;
    n++;
    y=y0+(2*j+1)*dy;
    for(i=0;i<x_size-1;i++) {
      mesh.vertices[n].lon=x0+i*dx+dx/2.0;
      mesh.vertices[n].lat=y;
      n++;
      }
    mesh.vertices[n].lon=x0;
    mesh.vertices[n].lat=y;
    n++;
    mesh.vertices[n].lon=x0+length;
    mesh.vertices[n].lat=y;
    n++;
    }
  
  mesh.nvtxs=n;
  
  for(n=0;n<mesh.nvtxs;n++) {
    for(m=0;m<mesh.nvtxs;m++) {
      if(n==m) continue;
      x=mesh.vertices[m].lon-mesh.vertices[n].lon;
      y=mesh.vertices[m].lat-mesh.vertices[n].lat;
      d=sqrt(x*x+y*y);
      if(d<dx*1.0001) {
        status=fe_connectvertices(mesh, m,n);
        }
      }
    }
/* *-----------------------------------------------------------------------------
  finalize mesh*/
  for (n=0;n<mesh.nvtxs;n++) {
    x=mesh.vertices[n].lon;
    y=mesh.vertices[n].lat;
//     projection_to_geo(sgrid.projection,&(mesh.vertices[n].lat),&(mesh.vertices[n].lon),x,y);
    status=geo_mercator_inverse(projection,&mesh.vertices[n].lon,&mesh.vertices[n].lat,x,y);   
    }
  status=fe_savemesh("hexagonal.nei",MESH_FILE_FORMAT_TRIGRID,mesh);
  return(0);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_QtoT(mesh_t & quadmesh, mesh_t & trimesh, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,m,n,status;
  double x,y;
  
  trimesh.nvtxs=quadmesh.nvtxs+quadmesh.nquadrangles;
  trimesh.ntriangles=4*quadmesh.nquadrangles;

  trimesh.vertices =new vertex_t[trimesh.nvtxs];
  trimesh.triangles=new triangle_t[trimesh.ntriangles];

/*------------------------------------------------------------------------------
  recycle quadrangle vertices */
  for(n=0;n<quadmesh.nvtxs;n++) {
    x=quadmesh.vertices[n].lon;
    y=quadmesh.vertices[n].lat;
    trimesh.vertices[n].lon=x;
    trimesh.vertices[n].lat=y;
/*------------------------------------------------------------------------------
    build connections */
    trimesh.vertices[n].nngh=quadmesh.vertices[n].nngh;
    trimesh.vertices[n].ngh=new int[trimesh.vertices[n].nngh];
    for(i=0;i<quadmesh.vertices[n].nngh;i++) {
      int nn=quadmesh.vertices[n].ngh[i];
      trimesh.vertices[n].ngh[i]=nn;
      }
    }
    
/*------------------------------------------------------------------------------
  and add quadrangle barycenter */
  n=quadmesh.nvtxs;
  for(m=0;m<quadmesh.nquadrangles;m++) {
    status=fe_position(quadmesh, quadmesh.quadrangles[m], &x, &y,0);
    trimesh.vertices[n].lon=x;
    trimesh.vertices[n].lat=y;
/*------------------------------------------------------------------------------
    build connections */
    trimesh.vertices[n].nngh=4;
    trimesh.vertices[n].ngh=new int[4];
    for(i=0;i<4;i++) {
      int nn=quadmesh.quadrangles[m].vertex[i];
      trimesh.vertices[n].ngh[i]=nn;
      }
    status=fe_insertvertex(&trimesh, n);
    n++;
    }
  
// /*------------------------------------------------------------------------------
//   build connections */
//   for(m=0;m<quadmesh.nquadrangles;m++) {
//     for(i=0;i<4;i++) {
//       n=quadmesh.quadrangles->vertex[i];
//       }
//     }
    
  status=fe_cleanvertices(&trimesh, false);
  
  status=fe_list(&trimesh);
  
  status=fe_edgetable(&trimesh,0,0);
  status=fe_codetable2(&trimesh,0,1,0);
    
  return(status);
  
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_create_rectangular(grid_t sgrid, grid_t cgrid, geo_t projection, criteria_t criteria, float *topo, float topomask, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  create various rectangular shape meshes:
  
  * fix resolution quadrangles
  * fix resolution triangles, square distribution
  * fix resolution triangles, hexagonal distribution
  * variable resolution triangles
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/  
{
  int i,j,k,m,status;
  char *landmask, *z_landmask;
  grid_t topogrid, z_grid;
  double width;
  int i1, i2, j1, j2;
  topo_t synthetic;
  double slope;
  mesh_t mesh, mimic, mimicT;
  plg_t *polygones=0;
  plg_t *limits=0;
  int nlimits;  
  char belname[1024];
  char *meshsize=NULL;
   
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
 UNIFORM RESOLUTION
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  create quadrangle mesh with notebook parameters (uniform resolution)*/
#if 0
  status=fe_create_rectangularQ(sgrid, sgrid, topo, topomask, debug);
#else
  int FixPinched=0;
  z_grid=map_f2zgrid(sgrid);
  z_landmask=new char [z_grid.Hsize()];
  for(size_t n=0;n<z_grid.Hsize();n++) {
//    z_landmask[n]=(ztopo[n]!=tag);
    z_landmask[n]=1;
    }
  status=fe_ZGrid2Quadrangle(z_grid, sgrid, z_landmask, topo, topomask, mimic, FixPinched, 0);
  status=fe_savemesh("uniformQ.nei",MESH_FILE_FORMAT_TRIGRID,mimic);
#endif
  
  status=fe_QtoT(mimic, mimicT, debug);
  status=fe_savemesh("uniformTbis.nei",MESH_FILE_FORMAT_TRIGRID,mimicT);
  
/*------------------------------------------------------------------------------
  create corresponding triangular mesh, square distribution */  
  z_grid=map_f2zgrid(sgrid);
  z_landmask=new char [z_grid.Hsize()];
  for(size_t n=0;n<z_grid.Hsize();n++) {
//    z_landmask[n]=(ztopo[n]!=tag);
    z_landmask[n]=1;
    }
  landmask=new char [sgrid.Hsize()];
  for(size_t n=0;n<sgrid.Hsize();n++) {
    landmask[n]=1;
    }
    
  status=fe_mimicgrid(sgrid, &mimic, landmask, topo, z_grid, z_landmask, (const char *) 0, debug); 
  status=fe_savemesh("uniformT.nei", MESH_FILE_FORMAT_TRIGRID, mimic);
  
  mimic.destroy();
  z_grid.free();
  
  delete[] z_landmask;
  delete[] landmask;
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
 VARIABLE RESOLUTION
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

/*------------------------------------------------------------------------------
  create quadrangle mesh with variable resolution */

  float *density=new float[cgrid.Hsize()], mask=1.e+10;
  
  for(size_t m=0; m<cgrid.Hsize(); m++)  density[m]=0.0;
  
//   criteria.minsize=0.5;
//   criteria.maxsize=0.5;
// 
//   criteria.shelf_minsize=0.5;
//   criteria.shelf_maxsize=0.5;
//   
//   criteria.mode=1+2+0*4+0*8+0*16+1*32;
//   
//   criteria.cellsize=0.05;
//   status=fe_save_criteria("dum.crt", criteria);

  status=fe_ComputeMeshsize("density.nc", criteria, &cgrid, &density, &mask, 0, 0, 0, false);
  
  double ymean=0.5*(sgrid.y[0]+sgrid.y[sgrid.nx]);
  plg_t target(sgrid.x[0],ymean,sgrid.x[sgrid.nx-1],ymean), sampled;
  
  criteria.cellsize=0.025;
  target.SetFlag('M');      
  
  status=plg_cartesian(cgrid.projection, &target, 1);

  sampled=plg_sample(target, cgrid, density, mask, criteria, 0.95, 0, false);
  
  status=plg_save("mesh-academic.variable.plg", PLG_FORMAT_SCAN,  &sampled, 1);

  grid_t variable;
  
  variable.init(sampled.npt,2);
  
  for(j=0; j<variable.ny;j++) {
    double p=sgrid.y[j*sgrid.nx+j];
    for(i=0; i<variable.nx;i++) {
      size_t m=variable.nx*j+i;
      variable.x[m]=sampled.t[i];
      variable.y[m]=p;
      }
    }
  
  z_grid=map_f2zgrid(variable);
  z_landmask=new char [z_grid.Hsize()];
  for(size_t n=0;n<z_grid.Hsize();n++) {
//    z_landmask[n]=(ztopo[n]!=tag);
    z_landmask[n]=1;
    }
  landmask=new char [variable.Hsize()];
  for(size_t n=0;n<variable.Hsize();n++) {
    landmask[n]=1;
    }
    
  status=fe_ZGrid2Quadrangle(z_grid, variable, z_landmask, topo, topomask, mimic, FixPinched, 0);
  status=fe_savemesh("variableQ.nei",MESH_FILE_FORMAT_TRIGRID,mimic);
  
  status=fe_mimicgrid(variable, &mimic, landmask, topo, z_grid, z_landmask, (const char *) 0, debug); 
  status=fe_savemesh("variableT.nei",MESH_FILE_FORMAT_TRIGRID,mimic);
  
  mimic.destroy();
  z_grid.free();
  
  delete[] z_landmask;
  delete[] landmask;
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
 HEXAGONAL SHAPE, UNIFORM RESOLUTION
 
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/*------------------------------------------------------------------------------
  create triangular mesh, hexagonal distribution */

  double dx=MIN(cgrid.dx,cgrid.dy);
  
  i=1;
  j=0;
  m=cgrid.nx*j+i;
  double x0=cgrid.x[m];
  double y0=cgrid.y[m];
  
  i=cgrid.nx-2;
  j=0;
  m=cgrid.nx*j+i;
  double length=cgrid.x[m]-x0;
  
  i=1;
  j=cgrid.ny-2;
  m=cgrid.nx*j+i;
  width=cgrid.y[m]-y0;
  
  width=10000.;
  
  int x_size=1+length/dx;
  int y_size=width/dx/2.0;
  y_size=2;
  dx=length/(x_size-1);

  status=fe_create_hexagonal(sgrid, cgrid, projection, criteria, dx, length, y_size, debug);
  
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
 
  UNDISCIPLINED SHAPE, VARIABLE RESOLUTION
 
XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

/*------------------------------------------------------------------------------
  create variable resolution mesh */
  nlimits=1;
    
  width=min(10000.,cgrid.dy);
  
  j1=cgrid.ny/2-floor(width/2.0/cgrid.dy+0.5);
  j2=cgrid.ny/2+floor(width/2.0/cgrid.dy+0.5);
  
//   i1=1;
//   i2=sgrid.nx-2;
  
  i1=0;
  i2=sgrid.nx-1;
  
  j1=max(0,j1);
  j2=min(cgrid.ny-1,j2);
  
  limits=new plg_t[nlimits];
  limits[0].npt=5;
  limits[0].t=new double[limits[0].npt];
  limits[0].p=new double[limits[0].npt];
  
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

  trouble with plg_cartesian_core_template called from 
  int plg_cartesian(projPJ ref, vector<plg_t> & polygones) if not initialized

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  limits[0].x=new double[limits[0].npt];
  limits[0].y=new double[limits[0].npt];
  
  k=0;
  i=i1;
  j=j1;
  m=sgrid.nx*j+i;
  limits[0].t[k]=sgrid.x[m];
  limits[0].p[k]=sgrid.y[m];
  
  k++;
  i=i2;
  j=j1;
  m=sgrid.nx*j+i;
  limits[0].t[k]=sgrid.x[m];
  limits[0].p[k]=sgrid.y[m];
  
  k++;
  i=i2;
  j=j2;
  m=sgrid.nx*j+i;
  limits[0].t[k]=sgrid.x[m];
  limits[0].p[k]=sgrid.y[m];
  
  k++;
  i=i1;
  j=j2;
  m=sgrid.nx*j+i;
  limits[0].t[k]=sgrid.x[m];
  limits[0].p[k]=sgrid.y[m];
  
  k++;
  i=i1;
  j=j1;
  m=sgrid.nx*j+i;
  limits[0].t[k]=sgrid.x[m];
  limits[0].p[k]=sgrid.y[m];
  
  if(debug) status=plg_save("mesh-academic.plg", PLG_FORMAT_SCAN,  limits, nlimits);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create synthetic bel file 

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  FILE  *fp;
  sprintf(belname,"mesh-academic.bel");
  printf("#################################################################\n");
  printf("BEL file generation: %s\n", belname);
  if ((fp = fopen (belname, "w")) == NULL ) {
    fprintf (stderr, "could not open file %s\n", belname);
    return(-1);
    }
  fprintf(fp,"%s \n","#OPEN BOUNDARIES");
  for (k=0; k<limits[0].npt-1;k++){
    if ((k==1) || (k==3)) {
//       fprintf(fp,"%lf %lf %lf %lf 5\n" ,  limits[0].t[k],limits[0].p[k],limits[0].t[k+1],limits[0].p[k+1]);
      fprintf(fp,"%g %g %g %g 5\n" ,  limits[0].t[k],limits[0].p[k],limits[0].t[k+1],limits[0].p[k+1]);
      }
    if ((k==0) || (k==2)) {
//       fprintf(fp,"%lf %lf %lf %lf 20\n" , limits[0].t[k],limits[0].p[k],limits[0].t[k+1],limits[0].p[k+1]);
      fprintf(fp,"%g %g %g %g 20\n" , limits[0].t[k],limits[0].p[k],limits[0].t[k+1],limits[0].p[k+1]);
      }
    }
  fclose (fp);
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create variable resolution mesh

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_create_rectangularT(criteria, mesh, limits, nlimits, meshsize, debug);

  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_create_COMODO(grid_t sgrid, grid_t cgrid, geo_t projection, criteria_t criteria, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  create various rectangular shape meshes:
  
  * fix resolution quadrangles
  * fix resolution triangles, square distribution
  * fix resolution triangles, hexagonal distribution
  * variable resolution triangles
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/  
{
  int i,j,k,m,status;
  float *topo, topomask=1.e+35;
  char *landmask, *z_landmask;
  grid_t topogrid, z_grid;
  double width;
  int j1, j2;
  topo_t synthetic;
  double slope;
  mesh_t mesh, mimic, mimicT;
  plg_t *polygones=0;
  plg_t *limits=0;
  int nlimits;  
  char belname[1024];
  char *meshsize=NULL;
  
/* Original setup */
/*------------------------------------------------------------------------------
                                                   length, slope, depth (m);  */
//   synthetic.shore.set     (TOPO_FIXED,        15.0e+03, 0.0,       10.0);
//   synthetic.shorebreak.set(TOPO_TRANSITION,    5.0e+03, 1.e-02,    10.0);
//   synthetic.shelf.set     (TOPO_FIXED,        75.0e+03, 0.0,      100.0);
//   synthetic.shelfbreak.set(TOPO_TRANSITION,   25.0e+03, 2.e-01,   100.0);
//   synthetic.deep.set      (TOPO_RESIZABLE,   500.0e+03, 0.0,     1000.0);
   
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  In COMODO IW test case slope must shelf breack from 200m to 4000m depth in 50km
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
  slope=3.8/50;
/*------------------------------------------------------------------------------
                                                length, slope, depth (m);  */
  synthetic.shore.set     (TOPO_FIXED,        15.0e+03,   0.0,    200.0);
  synthetic.shorebreak.set(TOPO_TRANSITION,    5.0e+03,   0.0,    200.0);
  synthetic.shelf.set     (TOPO_FIXED,       410.0e+03,   0.0,    200.0);
  synthetic.shelfbreak.set(TOPO_TRANSITION,   50.0e+03, slope,    200.0);
  synthetic.deep.set      (TOPO_RESIZABLE,   400.0e+03,   0.0,   4000.0);

  topo=new float[cgrid.nx*cgrid.ny];
  status=topo_synthetic(synthetic, "depth.nc", cgrid, sgrid, topo, topomask, false);
  
#if 1
  status=comodo_depth(cgrid, -200., -4000.0, topo, topomask);
  status=topo_savefield("comodo-depth.nc", sgrid, topo, topomask);
  float *smoothed=new float[sgrid.nx*sgrid.ny];
  float scale=20.0e+03;
  sgrid.dx=sgrid.x[1]-sgrid.x[0];
  sgrid.dy=sgrid.y[sgrid.nx]-sgrid.y[0];
  status=map_smooth_latThenLong(sgrid, topo, topomask, scale, smoothed);
  status=topo_savefield("comodo-smoothed.nc", sgrid, smoothed, topomask);
#endif
  
  status=fe_create_rectangular(sgrid, cgrid, projection, criteria, topo, topomask, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_create_GEOVIDE(grid_t sgrid, grid_t cgrid, geo_t projection, criteria_t criteria, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  create various rectangular shape meshes:
  
  * fix resolution quadrangles
  * fix resolution triangles, square distribution
  * fix resolution triangles, hexagonal distribution
  * variable resolution triangles
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/  
{
  int status;
  float *topo, topomask=1.e+35;
  grid_t topogrid, z_grid;
  shape_t synthetic;
  atom_shape_t *section;
  double slope;
  char belname[1024];
  string filename="topo-section-01.txt";
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  In COMODO IW test case slope must shelf breack from 200m to 4000m depth in 50km
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  slope=3.8/50;
/*------------------------------------------------------------------------------
                                                length, slope, depth (m);  */
//   synthetic.shelf.set     (TOPO_FIXED,         5.0e+03,   0.0,      2.0);
//   synthetic.shelfbreak.set(TOPO_TRANSITION,   50.0e+03, slope,    200.0);
//   synthetic.deep.set      (TOPO_RESIZABLE,   800.0e+03,   0.0,   3500.0);

/*------------------------------------------------------------------------------
                                          length, slope, depth (m) depth (m)  */
  section=new atom_shape_t;
  section->set_linear(TOPO_FIXED,          5000.0,   NAN,     10.0,   10.0);
  synthetic.transects.push_back(*section);
  section=new atom_shape_t;
  section->set_external(TOPO_SECTION_FILE, filename);
  synthetic.transects.push_back(*section);
  section=new atom_shape_t;
  section->set_linear(TOPO_RESIZABLE,   800.0e+03,   0.0,   3500.0,   3500.0);
  synthetic.transects.push_back(*section);

  topo=new float[cgrid.nx*cgrid.ny];
  status=topo_synthetic(synthetic, "depth.nc", cgrid, sgrid, topo, topomask, true);
    
  status=fe_create_rectangular(sgrid, cgrid, projection, criteria, topo, topomask, debug);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int fe_create_WaveCanal(grid_t sgrid, grid_t cgrid, geo_t projection, criteria_t criteria, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  create various rectangular shape meshes:
  
  * fix resolution quadrangles
  * fix resolution triangles, square distribution
  * fix resolution triangles, hexagonal distribution
  * variable resolution triangles
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/  
{
  int i,j,k,m,status;
  float *topo, topomask=1.e+35;
  char *landmask, *z_landmask;
  grid_t topo_cgrid, topo_sgrid, z_grid;
  double width;
  int j1, j2;
  topo_t synthetic;
  double slope, length;
  mesh_t mesh, mimic, mimicT;
  plg_t *polygones=0;
  plg_t *limits=0; 
  int nlimits;  
  char belname[1024];
  char *meshsize=NULL;
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  flat canal
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  length=cgrid.x[cgrid.Hsize()-1]-cgrid.x[0]*2;

/*------------------------------------------------------------------------------
                                                   length, slope, depth (m);  */
//   synthetic.beach.set     (TOPO_TRANSITION,           200,   0.0,    -10.0);
  synthetic.shore.set     (TOPO_FIXED,             length,   0.0,     20.0);
  
  topo=new float[cgrid.nx*cgrid.ny];
  status=topo_synthetic(synthetic, "depth.nc", cgrid, sgrid, topo, topomask, false);
  synthetic.destroy();
  
  synthetic.shore.set      (TOPO_FIXED,                500,   0.0,    -10.0);
  synthetic.shorebreak.set (TOPO_TRANSITION,          1000,   0.0,      0.0);
  synthetic.shelf.set      (TOPO_FIXED,             length,   0.0,     20.0);
  
  topo=new float[cgrid.nx*cgrid.ny];
  status=topo_synthetic(synthetic, "beach.nc", cgrid, sgrid, topo, topomask, false);
  
  status=fe_create_rectangular(sgrid, cgrid, projection, criteria, topo, topomask, debug);
  
  resize_t resize;
  resize.init(cgrid, -5, -5, -10, -10);
  status=resize_grid(cgrid, topo_cgrid, resize);
  topo_sgrid=map_get_spherical(topo_cgrid.projection,topo_cgrid);
  map_printgrid(sgrid);
  map_printgrid(topo_sgrid);
  
  delete[] topo;

  topo=new float[topo_cgrid.nx*topo_cgrid.ny];
  status=topo_synthetic(synthetic, "depth.nc", topo_cgrid, topo_sgrid, topo, topomask, false);
  
  return(status);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int mercator(grid_t & cgrid, grid_t & sgrid, geo_t projection)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status=0;
  char parameters[1024];
  projPJ pj=assign_projection(0, &projection.lat, &projection.lon, parameters);
  cgrid.projection=pj;
  cgrid.proj4options=poc_strdup(parameters);
  
  sgrid=map_get_spherical(pj,cgrid);
//   for(int m=0; m< cgrid.Hsize(); m++) {
//     cgrid.x[m]=sgrid.x[m];
//     cgrid.y[m]=sgrid.y[m];
//     }
//   status=geo_to_projection (cgrid.proj4options, cgrid.x, cgrid.y, cgrid.Hsize());
//   
//   status=map_minmax(&cgrid);  
  
  return(status);
}
/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int i,j,k,m;
  int option,channels=0;
  FILE *file,*in;
  FILE *out,*echo;
  char *keyword,*zone=NULL,*gridfile=NULL,*criteriafile=NULL;
  char *input=NULL,*meshfile=NULL,*depthfile=NULL,*output=NULL,*poly=NULL,*meshsize=NULL;
  char *boundaryfile=NULL,*descriptor=NULL;
  char *notebook=NULL;
  char neighfile[1024],additional[1024],belname[1024];
  criteria_t criteria;
  plg_array_t set;
  int format=PLG_FORMAT_BOUNDARIES,npolygones=0;
  float *density,mask;
  double x,y;
  grid_t sgrid,cgrid;
  geo_t projection;
  bool debug=false;
  
  fct_echo( argc, argv);

  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        if(strcmp(keyword,"--debug")==0) {
          debug=true;
          n++;
          continue;
          }
        switch (keyword[1]) {
        case 'b' :
          depthfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'c' :
          criteriafile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'd' :
          descriptor= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          boundaryfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'n' :
          notebook= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'p' :
          poly= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 's' :
          meshsize= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'z' :
          zone= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        input= strdup(argv[n]);
        n++;
        break;
      }
    free(keyword);
    }
    
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create north-pole radial mesh
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   status= fe_create_circularT(debug);
//   exit(0);
  
/*------------------------------------------------------------------------------
  load and initialize criteria structure*/
  status=fe_save_criteria("empty.criteria", criteria);
  if(criteriafile != NULL) {
    status=fe_load_criteria(criteriafile, &criteria);
    }

  ntbk_grid_t ngrid;
  status=save_notebook("empty.notebook", ngrid);

 if(notebook != NULL) {
    status=load_notebook(notebook, &cgrid, &sgrid, &projection);
    if(status !=0) {
      printf("cannot properly load notebook file=%s\n",notebook);
      goto error;
      }
    status=mercator(cgrid, sgrid, projection);
    printf("%s (notebook file) processed\n",notebook);
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create COMODO IW testcase mesh
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   status=fe_create_COMODO(sgrid, cgrid, projection, criteria, debug);
//   exit(0);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create non-hydrostatic wave canal mesh
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   status=fe_create_WaveCanal(sgrid, cgrid, projection, criteria, debug);
//   exit(0);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create SHOM 2018 academic estuary
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

//   status=fe_create_Estuary(debug);
//   exit(0);

/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  create GEOVIDE 2019 internal tide mesh
 
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  status=fe_create_GEOVIDE(sgrid, cgrid, projection, criteria, debug);
  exit(0);

  if(status !=0) {
    printf("mesh generation failed...\n");
    goto error;
    }
  
  printf("academic mesh sucessfully completed\n");

end: __OUT_BASE_LINE__("end of mesh-academic ... \n");
  exit(0);

error:
  __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}
