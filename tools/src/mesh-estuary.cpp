

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

/**-************************************************************************

  2D mesh generator (academic, simple geometry meshes)

  F. Lyard, 2009, CNRS/LEGOS, Toulouse, France

  e-mail: florent.lyard@legos.obs-mip.fr

***************************************************************************/

// #include <stdio.h>
// #include <stdlib.h>
// #include <string.h>
// #include <list>
// 
// using namespace std;
// 
// #include "tools-structures.h"
// 
// #include "functions.h"
// #include "geo.h"
// #include "fe.h"
// #include "archive.h"
// #include "map.h"
// #include "polygones.h"
// #include "grd.h"
// #include "list.h"
// #include "maths.h"
// #include "sym-io.h"
// #include "map.def"
// #include "mass.h"
// #include "topo.h"
// 
// #include "zapper.h"     /*  rutin.h contains common utility routines  */
// 
// extern  void fe_LGP0_to_LGP1(double *in, double *out, mesh_t mesh);
// #include "statistic.h"
// 
// #define TOPO_UNDEFINED   -1
// #define TOPO_FIXED        0
// #define TOPO_TRANSITION   1
// #define TOPO_RESIZABLE    2
// 
// class shape_t {
// private :
// public :
//   int type;
// 
//   shape_t() {
//     type=-1;
//     }
// 
//   void reset() {
//     type=-1;
//     }
// 
//   void destroy() {
//     (*this).reset();
//     }
// };
// 
// class topo_section_t {
// private :
// public :
//   int type;
//   double length, slope, h0;
//   double limits[2];
// 
//   topo_section_t() {
//     type=TOPO_UNDEFINED;
//     length=slope=h0=0;
//     limits[0]=limits[1]=-1.;
//     }
// 
//   void set(int t, double l, double a, double h) {
//     type=t;
//     length=l;
//     slope=a;
//     h0=h;
//     }
// 
//   void reset() {
//     type=-1;
//     }
// 
//   void destroy() {
//     (*this).reset();
//     }
//     
//   int depth(double d, double & depth) {
//     double x=0;
//     int status;
//     
//     d=d-limits[0];
//     
//     if(length>0) {
//       x=d/length;
//       if( (x>= 0.0) && (x<=1.0)) {
//         depth=h0+slope*d;
//         status=0;
//         }
//       else {
//         status=-1;
//         }
//       }
//     else {
//       status=-1;
//       }
//     return(status);
//     }
// 
// };
// 
// class topo_t {
// private :
// public :
//   int proportional,type;
//   topo_section_t ridge,deep,shelfbreak,shelf,shorebreak,shore;
// 
//   topo_t() {
//     type=-1;
//     proportional=0;
//     }
// 
//   void reset() {
//     type=-1;
//     }
// 
//   void destroy() {
//     (*this).reset();
//     }
//     
//   int depth(double d, double & depth) {
//     double L=0;
//     int status;
//     list<topo_section_t>::iterator it;
//     
//     list<topo_section_t> list;
//     list.push_back(shore);
//     list.push_back(shorebreak);
//     list.push_back(shelf);
//     list.push_back(shelfbreak);
//     list.push_back(deep);
//     list.push_back(ridge);
//     
//     it=list.begin();
//     do {
//       status=(*it).depth(d,depth);
//       it++;
//       } while( (status==-1) && (it != list.end()) );
//     type=-1;
//     list.clear();
//     
//     return(status);
//     }
// 
//   double length(double & resizable) {
//     double d, L=0;
//     double depth;
//     int status;
//     vector<topo_section_t>::iterator it,previous;
//     
//     vector<topo_section_t> list;
//     list.push_back(shore);
//     list.push_back(shorebreak);
//     list.push_back(shelf);
//     list.push_back(shelfbreak);
//     list.push_back(deep);
//     list.push_back(ridge);
//     
//     resizable=0;
//     
//     it=list.begin();
//     do {
//       L+=(*it).length;
//       if((*it).type==TOPO_RESIZABLE) resizable+=(*it).length;
//       it++;
//       } while( it != list.end() );
//       
//     shore=list[0];
//     shorebreak=list[1];
//     shelf=list[2];
//     shelfbreak=list[3];
//     deep=list[4];
//     
//     list.clear();
//     return(L);
//     }
// 
//   void resize(double scale) {
//     double d, L=0;
//     double depth;
//     int status;
//     vector<topo_section_t>::iterator it,previous;
//     
//     vector<topo_section_t> list;
//     list.push_back(shore);
//     list.push_back(shorebreak);
//     list.push_back(shelf);
//     list.push_back(shelfbreak);
//     list.push_back(deep);
//     list.push_back(ridge);
//         
//     it=list.begin();
//     do {
//       if((*it).type==TOPO_RESIZABLE) (*it).length*=scale;
//       it++;
//       } while( it != list.end() );
//       
//     shore=list[0];
//     shorebreak=list[1];
//     shelf=list[2];
//     shelfbreak=list[3];
//     deep=list[4];
//     
//     list.clear();
//     }
// 
//   int parse(double size) {
//     double d, L=0, resizable, scale;
//     double depth;
//     int status;
//     vector<topo_section_t>::iterator it,previous,next;
//     
//     vector<topo_section_t> list;
//     
// /**----------------------------------------------------------------------------
//     adjust lengths */
//     L=this->length(resizable);
//     
//     scale=(size+resizable-L)/resizable;
//     
//     this->resize(scale);
//     
//     L=this->length(resizable);
// 
//     list.push_back(shore);
//     list.push_back(shorebreak);
//     list.push_back(shelf);
//     list.push_back(shelfbreak);
//     list.push_back(deep);
//     list.push_back(ridge);
// 
//     it=list.begin();
//     status=(*it).limits[0]=0.0;
//     status=(*it).limits[1]=(*it).length;
//     do {
//       previous=it;
//       it++;
//       status=(*it).limits[0]=(*previous).limits[1];
//       status=(*it).limits[1]=(*previous).limits[1]+(*it).length;
//       } while( it != list.end() );
//     
// /**----------------------------------------------------------------------------
//     adjust depths */
//     it=list.begin();
//     it++;
//     do {
//       if((*it).type==TOPO_TRANSITION) {
//         previous=next=it;
//         previous--;
//         next++;
//         double h[2];
//         status=(*it).depth((*it).limits[0],h[0]);
//         status=(*it).depth((*it).limits[1],h[1]);
//         status=(*previous).depth((*it).limits[0],h[0]);
//         status=(*next).depth    ((*it).limits[1],h[1]);
//         (*it).h0=h[0];
//         (*it).slope=(h[1]-h[0])/(*it).length;
//         }
//       it++;
//       } while( it != list.end() );
//     
// 
//     shore=list[0];
//     shorebreak=list[1];
//     shelf=list[2];
//     shelfbreak=list[3];
//     deep=list[4];
//     
//     list.clear();
//     return(0);
//     }
// };
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// template <typename TYPE> int topo_synthetic(topo_t synthetic, const char *bathymetry, grid_t grid, grid_t sgrid, TYPE *topo, TYPE topomask)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int   i,j,k,l,m,n,status;
//   int   ii,jj,imin,imax,jmin,jmax;
//   float *buffer,*radius,mask;
//   double z,d, maxscale;
//   signed char   *selected;
//   double x,y,dx,dy;
//   
//   synthetic.parse(grid.x[grid.nx-1]-grid.x[0]);
// 
//   for(j=0;j<grid.ny;j++) {
//     for(i=0;i<grid.nx;i++) {
//       m=j*grid.nx+i;
//       topo[m]=topomask;
//       }
//     }
// 
// /* *------------------------------------------------------------------------------
//   */
//   for(j=0;j<grid.ny;j++) {
//     for(i=0;i<grid.nx;i++) {
//       m=j*grid.nx+i;
//       d=grid.x[m]-grid.x[0];
//       status=synthetic.depth(d,z);
//       if(status==0) {
//         topo[m]=-z;
//         }
//       else {
//         status=synthetic.depth(d,z);
//         }
//       }
//     }
//     
//   status=topo_savefield(bathymetry, sgrid, topo, topomask);
//  
//   TYPE *smoothed=new TYPE[sgrid.nx*sgrid.ny];
//   float scale=20.0e+03;
//   
//   sgrid.dx=sgrid.x[1]-sgrid.x[0];
//   sgrid.dy=sgrid.y[sgrid.nx]-sgrid.y[0];
//   status=map_smooth_latThenLong(sgrid, topo, topomask, scale, smoothed);
// 
//   status=topo_savefield("smoothed.nc", sgrid, smoothed, topomask);
//  
//   return(0);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int trisolverd(int la, double *a, double *b, double *diag, double *rhs)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//     /* Local variables */
//     int k;
//     double tmp, tmpdiag;
// 
// /* -------------------------------------------------------------- */
// /* TRIAGONAL SYSTEM SOLVER */
// /* -------------------------------------------------------------- */
// /* first Gauss elimination */
// 
// /* At surface */
// 
//     /* Function Body */
//     tmp = 1. / diag[0];
//     a[0] *= tmp;
//     rhs[0] *= tmp;
// 
// /* at mid depth */
// 
//     for (k = 1; k < la-1; k++) {
//         tmpdiag = diag[k] - a[k - 1] * b[k];
//         tmp = 1. / tmpdiag;
//         a[k] *= tmp;
//         rhs[k] = (rhs[k] - rhs[k - 1] * b[k]) * tmp;
//     }
// 
// /* at bottom, solution in rhs */
// 
//     tmpdiag = diag[la-1] - a[la-2] * b[la-1];
//     tmp = 1. / tmpdiag;
//     rhs[la-1] = (rhs[la-1] - rhs[la-2] * b[la-1]) * tmp;
// /* -------------------------------------------------------------- */
// /* second and final Gauss elimination to surface */
// /* solution in rhs */
// 
//     for (k = la-2; k >= 0; k--) {
//         rhs[k] -= rhs[k + 1] * a[k];
//     }
//     return 0;
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int integration(int n, double dx, double *rhs)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status; 
//   double *a, *b, *diag;
//   double dx2=dx*dx;
//   
//   a=new double[n];
//   b=new double[n];
//   
//   diag=new double[n];
//   for(int i=0;i<n;i++) {
//     a[i]=0;
//     b[i]=0;
//     diag[i]=0;
//     }
// 
//   diag[0]  =1.0;
//   diag[n-1]=1.0;
// 
//   for(int i=1;i<n-1;i++) {
//     a[i]=1.0/dx2;
//     b[i]=1.0/dx2;
//     diag[i]=-2.0/dx2;
//     }
//   
//   status=trisolverd(n, a, b, diag, rhs);
//   
//   delete[] a;
//   delete[] b;
//   delete[] diag;
// 
//   return 0;
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int comodo_depth(grid_t cgrid, float h0, float h1, float *topo, float mask)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status; 
//   int m, n;
//   double x, dx,  *rhs;
//   FILE *out;
//   
//   n=cgrid.nx;
//   dx=cgrid.dx;
//   
//   rhs=new double[n];
//   
//   for(int i=0; i<n; i++) {
//     m=i;
//     x=cgrid.x[m];
//     rhs[i]=1.e-08;
//     }
//     
//   rhs[0]  =h0;
//   rhs[n-1]=h1;
//   
//   status=integration(n, dx, rhs);
//   
//   for(int j=0; j<cgrid.ny; j++) {
//     for(int i=0; i<cgrid.nx; i++) {
//       m=cgrid.nx*j+i;
//       topo[m]=rhs[i];
//       }
//     }
//     
//   out=fopen("comodo-depth.gnu","w");
//   for(int i=0; i<n; i++) {
//     fprintf(out, "%lf %lf\n", cgrid.x[i], rhs[i]);
//     }
//   fclose(out);
//   
//   delete[] rhs;
// 
//   return 0;
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int fe_create_boundaries(plg_array_t *boundaries, plg_array_t plg, double radius)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status=0;
//   int k;
//   int equalize=0;
// 
//   boundaries->p=new plg_t[plg.n];
//   for(k=0;k<plg.n;k++) {
//     boundaries->p[k]=plg_resample(plg.p[k],radius, equalize);
//     }
//   return(status);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int fe_create_interiors(plg_array_t boundaries, double radius)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   __CHKERR_LINE_FILE__(ENOEXEC,"not finished");exit(ENOEXEC);
//   return(status);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int plg_create_rectangle(plg_array_t *plg, shape_t shape)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status=0;
//   int i,j,k;
//   double x[4],y[4];
//   double L,H,dt,dp,t0,p0,radius;
//   static double cycle[4][2]={(0,0),(1,0),(1,1),(0,1)};
// 
//   for(k=0;k<4;k++) {
//     x[k]=t0+cycle[k][0]*dt;
//     y[k]=p0+cycle[k][1]*dp;
//     }
// 
//   plg->n=4;
//   plg->p=new plg_t[plg->n];
// 
//   for(k=0;k<4;k++) {
//     i=k;
//     j=(k+1)%4;
//     plg->p[k]=plg_t(x[i],y[i],x[j],y[j]);
//     }
// 
//   return(status);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int plg_create_circle(plg_array_t *plg, shape_t shape, double radius, double t0, double p0)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int n,status;
//   int i,j,k,zone,ntoken;
//   double *x,*y,x0,y0;
//   double L,H;
//   char *proj4_options;
//   const char *hemisphere;
//   projPJ proj;
//   char **parms;
//   
//   L=2.*M_PI*radius;
// 
//   plg->n=1;
//   plg->p=new plg_t[plg->n];
// 
//   plg->p[0].npt=100;
//   
//   x=new double[plg->p[0].npt];
//   y=new double[plg->p[0].npt];
//   
//   if(p0<0) hemisphere="+south";
//   else hemisphere="";
//   
//   zone=(t0+180)/3;
//   
//   asprintf(&proj4_options,"+proj=utm +zone=%d +ellps=WGS84 +datum=WGS84 +units=m +no_defs %s",zone,hemisphere);
// 
//   proj = pj_init_plus(proj4_options);
//   if (!proj) __TRAP_ERR_EXIT__(1,"Projection initialization failed\n");
//   
//   geo_to_projection(proj, p0, t0, &x0, &y0);
//   pj_free(proj);
//   
//   for(n=0;n<plg->p[0].npt;n++) {
//     double alpha=2*M_PI/(plg->p[0].npt-1.);
//     x[n]=x0+radius*cos(alpha);
//     y[n]=x0+radius*sin(alpha);
//     }
//     
//   status=projection_to_geo (proj4_options, x, y, plg->p[0].npt);
//   
//   free(proj4_options);
// 
//   for(n=0;n<plg->p[0].npt;n++) {
//     plg->p[0].t[n]=x[n];
//     plg->p[0].t[n]=y[n];
//     }
// 
//   return(status);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// int plg_create_shape(plg_array_t *plg, shape_t shape)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   status=plg_create_rectangle(plg, shape);
// 
//   return(status);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int fe_create_rectangularT(criteria_t criteria, mesh_t & mesh, plg_t *limits, int nlimits, char *meshsize, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int i,j,k,l,m,n,status;
//   int count,s;
//   mesh_t *spherical;
//   double L,H,dt,dp,t0,p0,radius;
//   plg_array_t plg,boundaries;
//   shape_t shape;
//   point_t *points;
//   float *density,mask;
//   grid_t grid;
//   plg_t *polygones=0;
//   int npolygones=0;
//   const char *output="density.nc";
// //   bool debug=true;
//   
// // /* *----------------------------------------------------------------------------
// //   create a geometry (sort of synthetic "shorelines") */
// //   status=plg_create_shape(&plg, shape);
//   
// /**-----------------------------------------------------------------------------
//   echo criteria file */
//   status= fe_save_criteria("echo.crt", criteria);
//   status= fe_save_criteria("mesh-academic-echo.crt", criteria);
// 
// /* *-----------------------------------------------------------------------------
//   mesh size mapping*/
//   printf("#################################################################\n");
//   printf("mesh resolution setup\n");
//   if(meshsize==0) {
//     printf("compute mesh density\n");
//     status=fe_ComputeMeshsize(output, criteria, &grid, &density, &mask, limits, nlimits,1, false);
//     }
//   else {
//     printf("load mesh density from %s\n",meshsize);
//     status=fe_ReloadMeshsize(meshsize, &grid, &density, &mask, limits, nlimits);
//     status= defproj(&grid, limits, nlimits);
//     }
//     
//   printf("#################################################################\n");
//   printf("mesh limits sampling\n");
//   
//   for(int s=0; s<  nlimits; s++) {
//     }
//     
//   for(int s=0; s<  nlimits; s++) {
//     if(limits[s].npt==0) {
//       printf("empty polygon %d\n",s);
//       continue;
//       }
//     limits[s].flag=new char[limits[s].npt-1];
//     for(int k=0;k<limits[s].npt-1;k++) limits[s].flag[k]='M';
// 
//     plg_t *p=new plg_t[limits[s].npt-1];
//     plg_t *q=new plg_t[limits[s].npt-1];
//     for(int k=0;k<limits[s].npt-1;k++) {
//       p[k].duplicate(limits[s],k,k+1);
// /* *-----------------------------------------------------------------------------
//       necessary to insure mesh periodicity, to be improved*/
//       if(k==2) status=plg_flip(p[k]);
//       q[k]=plg_sample(p[k], grid, density, mask, criteria, 1, debug);
//       }
//     
// /* *-----------------------------------------------------------------------------
//     alternatively this should be done by applying a translation*/
// //     int N=q[2].npt-1;
// //     for(size_t l=1;l<N;l++) {
// //       q[2].t[l]=q[0].t[N-l]; //lon
// //       q[2].x[l]=q[0].x[N-l]; // x
// //       }
// 
//     for(int k=1;k<limits[s].npt-1;k++) plg_concat(q[0],q[k]);
// 
//     limits[s].duplicate(q[0]);
//     }
//   status=plg_save("mesh-academic.adjusted.plg", PLG_FORMAT_SCAN,  limits, nlimits);
// 
// // /* *----------------------------------------------------------------------------
// //   create mesh boundaries on the geometry */
// //   status=fe_create_boundaries(&boundaries, plg, radius);
// 
// /* *-----------------------------------------------------------------------------
//   mesh generation*/
//   printf("#################################################################\n");
//   printf("mesh generation\n");
//   mesh=fe_nodit(criteria, grid, density, mask, limits, nlimits, debug);
//   status= fe_savemesh("mesh-academic-no-reshape.nei",MESH_FILE_FORMAT_TRIGRID, mesh);
// 
//   return(0);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
// template <typename TYPE> int fe_create_rectangularQ(grid_t grid,grid_t topogrid, TYPE *topo, TYPE topomask)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
//   int status;
//   size_t count=0;
//   char *landmask;
//   mesh_t mesh;
//   double t,p;
//   
//   landmask=new char[grid.nx*grid.ny];
//  
//   for (size_t j=0;j<grid.ny;j++) {
//     for (size_t i=0;i<grid.nx; i++) {
//       float h;
//       double x=map_grid_x(grid,i,j);
//       double y=map_grid_y(grid,i,j);
//       x=map_recale(topogrid,x);
//       status=map_interpolation(topogrid,topo,topomask,x,y,&h);
//       size_t m=grid.nx*j+i;
//       if(h>0) {
//         landmask[m]=0;
//         }
//       else {
//         landmask[m]=1;
//         count++;
//         }
//       }
//     }
// 
//   status=fe_FGrid2Quadrangle(grid, landmask, topo,topomask,mesh,"quadrangle.nc",debug);
// 
//   status=fe_vertex_crosstables02(&mesh);
// 
//   int stopon_EdgeError =1;
//   int stopon_PinchError=1;
//   status=fe_codetable1(&mesh, 0, stopon_EdgeError, stopon_PinchError);
// 
//   status=fe_savemesh("test.nei",  MESH_FILE_FORMAT_TRIGRID, mesh);
//   
//   t=+5.;
//   p=-2.;
//   int n1=fe_nearest_vertex (mesh,t,p,0,0);
// 
//   t=+5.;
//   p=+2.;
//   int n2=fe_nearest_vertex (mesh,t,p,0,0);
//   
//   status=fe_setOBCflags(mesh, n1, n2, MESH_ELEVATION_NODE);
//   status=fe_write_boundarycode("test.bel", mesh, 1);
//   
//   return(0);
// }
// 
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int fe_create_hexagonal(grid_t sgrid, grid_t cgrid, geo_t projection, criteria_t criteria, double dx, double length, int y_size, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// /*------------------------------------------------------------------------------
//   
//   create rectangular shape, fix resolution, hexagonal distribution mesh
//     
// ------------------------------------------------------------------------------*/  
// {
//   int i,j,k,m,n,status;
//   mesh_t mesh;
//   double x0,y0,x,y,d;
//   
// //   double dx=criteria.shelf_minsize;
//   
//   i=1;
//   j=0;
//   m=cgrid.nx*j+i;
//   x0=cgrid.x[m];
//   
//   i=1;
//   j=cgrid.ny/2;
//   m=cgrid.nx*j+i;
//   y0=cgrid.y[m];
// //   
// //   i=cgrid.nx-2;;
// //   j=0;
// //   m=cgrid.nx*j+i;
// //   length=cgrid.x[m]-x0;
// //   
//   int x_size=1+length/dx;
//   
//   dx=length/(x_size-1);
//   
//   
//   double dy=dx*sqrt(3.0)/2.0;
//   
//   int half=y_size/2;
//   n=x_size+2*half*x_size;
//   
//   half=(y_size+1)/2;
//   n+=+2*half*(x_size+1);
//   
//   mesh.vertices=new vertex_t[n];
//   
//   n=0;
//   for(i=0;i<x_size;i++) {
//     mesh.vertices[n].lon=x0+i*dx;
//     mesh.vertices[n].lat=y0;
//     n++;
//     }
//     
//   half=y_size/2;  
//   for(j=0;j<half;j++) {
//     y=y0-2*(j+1)*dy;
//     for(i=0;i<x_size;i++) {
//       mesh.vertices[n].lon=x0+i*dx;
//       mesh.vertices[n].lat=y;
//       n++;
//       }
//     y=y0+2*(j+1)*dy;
//     for(i=0;i<x_size;i++) {
//       mesh.vertices[n].lon=x0+i*dx;
//       mesh.vertices[n].lat=y;
//       n++;
//       }
//     }
//   
//   half=(y_size+1)/2;  
//   for(j=0;j<half;j++) {
//     y=y0-(2*j+1)*dy;
//     for(i=0;i<x_size-1;i++) {
//       mesh.vertices[n].lon=x0+i*dx+dx/2.0;
//       mesh.vertices[n].lat=y;
//       n++;
//       }
//     mesh.vertices[n].lon=x0;
//     mesh.vertices[n].lat=y;
//     n++;
//     mesh.vertices[n].lon=x0+length;
//     mesh.vertices[n].lat=y;
//     n++;
//     y=y0+(2*j+1)*dy;
//     for(i=0;i<x_size-1;i++) {
//       mesh.vertices[n].lon=x0+i*dx+dx/2.0;
//       mesh.vertices[n].lat=y;
//       n++;
//       }
//     mesh.vertices[n].lon=x0;
//     mesh.vertices[n].lat=y;
//     n++;
//     mesh.vertices[n].lon=x0+length;
//     mesh.vertices[n].lat=y;
//     n++;
//     }
//   
//   mesh.nvtxs=n;
//   
//   for(n=0;n<mesh.nvtxs;n++) {
//     for(m=0;m<mesh.nvtxs;m++) {
//       if(n==m) continue;
//       x=mesh.vertices[m].lon-mesh.vertices[n].lon;
//       y=mesh.vertices[m].lat-mesh.vertices[n].lat;
//       d=sqrt(x*x+y*y);
//       if(d<dx*1.0001) {
//         status=fe_connectvertices(mesh, m,n);
//         }
//       }
//     }
// /* *-----------------------------------------------------------------------------
//   finalize mesh*/
//   for (n=0;n<mesh.nvtxs;n++) {
//     x=mesh.vertices[n].lon;
//     y=mesh.vertices[n].lat;
// //     projection_to_geo(sgrid.projection,&(mesh.vertices[n].lat),&(mesh.vertices[n].lon),x,y);
//     status=geo_mercator_inverse(projection,&mesh.vertices[n].lon,&mesh.vertices[n].lat,x,y);   
//     }
//   status=fe_savemesh("hexagonal.nei",MESH_FILE_FORMAT_TRIGRID,mesh);
//   return(0);
// }
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// 
//   int fe_create_rectangular(grid_t sgrid, grid_t cgrid, geo_t projection, criteria_t criteria, bool debug)
// 
// /*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
// {
// /*------------------------------------------------------------------------------
//   
//   create various rectangular shape meshes:
//   
//   * fix resolution quadrangles
//   * fix resolution triangles, square distribution
//   * fix resolution triangles, hexagonal distribution
//   * variable resolution triangles
//   
// ------------------------------------------------------------------------------*/  
//   int i,j,k,m,status;
//   float *topo, topomask=1.e+35;
//   char *landmask, *z_landmask;
//   grid_t topogrid, z_grid;
//   double width;
//   int j1, j2;
//   topo_t synthetic;
//   double slope;
//   mesh_t mesh, mimic;
//   plg_t *polygones=0;
//   plg_t *limits=0;
//   int nlimits;  
//   char belname[1024];
//   char *meshsize=NULL;
//   
// /* Original setup */
// //                                            Length , slope, depth in meters
// //   synthetic.shore.set     (TOPO_FIXED,        15.0e+03, 0.0,     10.0);
// //   synthetic.shorebreak.set(TOPO_TRANSITION,    5.0e+03, 1.e-02,  10.0);
// //   synthetic.shelf.set     (TOPO_FIXED,        75.0e+03, 0.0,    100.0);
// //   synthetic.shelfbreak.set(TOPO_TRANSITION,   25.0e+03, 2.e-01, 100.0);
// //   synthetic.deep.set      (TOPO_RESIZABLE,   500.0e+03, 0.0,   1000.0);
//    
//   
//   /************************
//    * In Comodo IW test case
//    * Slope must shelf breack from
//    * 200m to 4000m depth in 50km
//    */
//   slope=3.8/50;
//   synthetic.shore.set     (TOPO_FIXED,        15.0e+03, 0.0,    200.0);//Length , slope, depth in meters
//   synthetic.shorebreak.set(TOPO_TRANSITION,    5.0e+03, 0.0,    200.0);
//   synthetic.shelf.set     (TOPO_FIXED,       410.0e+03, 0.0,    200.0);
//   synthetic.shelfbreak.set(TOPO_TRANSITION,   50.0e+03, slope,  200.0);
//   synthetic.deep.set      (TOPO_RESIZABLE,   400.0e+03, 0.0,   4000.0);
// 
//   topo=new float[cgrid.nx*cgrid.ny];
//   status=topo_synthetic(synthetic, "depth.nc", cgrid, sgrid, topo, topomask);
// #if 1
//   status=comodo_depth(cgrid, 0.0, -4000.0, topo, topomask);
// #endif  
// /*------------------------------------------------------------------------------
//   create quadrangle mesh */
// #if 1
//   status=fe_create_rectangularQ(sgrid, sgrid, topo, topomask);
// #else
//   int FixPinched=0;
//   z_grid=map_f2zgrid(sgrid);
//   z_landmask=new char [z_grid.Hsize()];
//   for(size_t n=0;n<z_grid.Hsize();n++) {
// //    z_landmask[n]=(ztopo[n]!=tag);
//     z_landmask[n]=1;
//     }
//   status=fe_ZGrid2Quadrangle(z_grid, sgrid, z_landmask, topo, topomask, mimic, FixPinched, 0);
//   status=fe_savemesh("mimicQ.nei",MESH_FILE_FORMAT_TRIGRID,mimic);
// #endif
//   
// /*------------------------------------------------------------------------------
//   create triangular mesh, square distribution */  
//   z_grid=map_f2zgrid(sgrid);
//   z_landmask=new char [z_grid.Hsize()];
//   for(size_t n=0;n<z_grid.Hsize();n++) {
// //    z_landmask[n]=(ztopo[n]!=tag);
//     z_landmask[n]=1;
//     }
//   landmask=new char [sgrid.Hsize()];
//   for(size_t n=0;n<sgrid.Hsize();n++) {
//     landmask[n]=1;
//     }
//   
// //   debug=true;
//   
//   status=fe_mimicgrid(sgrid, &mimic, landmask, topo, z_grid, z_landmask, debug); 
//   status=fe_savemesh("mimicT.nei",MESH_FILE_FORMAT_TRIGRID,mimic);
// 
// /*------------------------------------------------------------------------------
//   create triangular mesh, hexagonal distribution */
// 
//   double dx=MIN(cgrid.dx,cgrid.dy);
//   
//   i=1;
//   j=0;
//   m=cgrid.nx*j+i;
//   double x0=cgrid.x[m];
//   double y0=cgrid.y[m];
//   
//   i=cgrid.nx-2;
//   j=0;
//   m=cgrid.nx*j+i;
//   double length=cgrid.x[m]-x0;
//   
//   i=1;
//   j=cgrid.ny-2;
//   m=cgrid.nx*j+i;
//   width=cgrid.y[m]-y0;
//   
//   width=10000.;
//   
//   int x_size=1+length/dx;
//   int y_size=width/dx/2.0;
//   
//   dx=length/(x_size-1);
// 
//   status=fe_create_hexagonal(sgrid, cgrid, projection, criteria, dx, length, y_size, debug);
//   
// /*------------------------------------------------------------------------------
//   create variable resolution mesh */
//   nlimits=1;
//     
//   width=10000.;
//   
//   j1=sgrid.ny/2-floor(width/2.0/cgrid.dy+0.5);
//   j2=sgrid.ny/2+floor(width/2.0/cgrid.dy+0.5);
//   
//   limits=new plg_t[nlimits];
//   limits[0].npt=5;
//   limits[0].t=new double[limits[0].npt];
//   limits[0].p=new double[limits[0].npt];
//   
//   k=0;
//   i=1;
//   j=j1;
//   m=sgrid.nx*j+i;
//   limits[0].t[k]=sgrid.x[m];
//   limits[0].p[k]=sgrid.y[m];
//   
//   k++;
//   i=sgrid.nx-2;
//   j=j1;
//   m=sgrid.nx*j+i;
//   limits[0].t[k]=sgrid.x[m];
//   limits[0].p[k]=sgrid.y[m];
//   
//   k++;
//   i=sgrid.nx-2;
//   j=j2;
//   m=sgrid.nx*j+i;
//   limits[0].t[k]=sgrid.x[m];
//   limits[0].p[k]=sgrid.y[m];
//   
//   k++;
//   i=1;
//   j=j2;
//   m=sgrid.nx*j+i;
//   limits[0].t[k]=sgrid.x[m];
//   limits[0].p[k]=sgrid.y[m];
//   
//   k++;
//   i=1;
//   j=j1;
//   m=sgrid.nx*j+i;
//   limits[0].t[k]=sgrid.x[m];
//   limits[0].p[k]=sgrid.y[m];
//   
//   if(debug) status=plg_save("mesh-academic.plg", PLG_FORMAT_SCAN,  limits, nlimits);
// 
// /*------------------------------------------------------------------------------
//   create synthetic bel file */  
//   FILE  *fp;
//   sprintf(belname,"mesh-academic.bel");
//   printf("#################################################################\n");
//   printf("BEL file generation: %s\n", belname);
//   if ((fp = fopen (belname, "w")) == NULL ) {
//     fprintf (stderr, "could not open file %s\n", belname);
//     return(-1);
//     }
//   fprintf(fp,"%s \n","#OPEN BOUNDARIES");
//   for (k=0; k<limits[0].npt-1;k++){
//     if ((k==1) || (k==3)) {
//       fprintf(fp,"%lf %lf %lf %lf 5\n" , limits[0].t[k],limits[0].p[k],limits[0].t[k+1],limits[0].p[k+1]);
//       }
//     }
//   fclose (fp);
//   
//   status=fe_create_rectangularT(criteria, mesh, limits, nlimits, meshsize, debug);
// 
//   return(status);
// }

