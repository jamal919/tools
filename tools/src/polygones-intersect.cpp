

#include <stdio.h>

#include "tools-structures.h"

#include "constants.h"
#include "polygones.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int plg_single(const plg_t & polygones, double lon, double lat, int *in, int mode, int check)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  M=A1 + l A1B1= A2 + l A2AB2

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int size,pos;
  int k,l;

  double *dx=NULL,*dy=NULL,*dn=NULL,*delta=NULL;
  double epsilon,error,tx,ty,ang,nx,ny,ux,uy,vx,vy;
  double alpha,beta;
  double *xp,*yp;
  bool colocated=false;
  
  switch(mode) {
    case PLG_SPHERICAL:
      xp=polygones.t;
      yp=polygones.p;
      break;
    case PLG_CARTESIAN:
      xp=polygones.x;
      yp=polygones.y;
      break;
    }
  
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx
 
  count the number of intersections of a half-line starting at (lon,lat) with an
  arbitrary direction and the polygone's segments
  
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/

  size=polygones.npt;

  exitIfNull(
    dx=new double[size]
    );
  exitIfNull(
    dy=new double[size]
    );
  exitIfNull(
    dn=new double[size]
    );
  exitIfNull(
    delta=new double[size]
    );
 
  error=1.e-05;
  epsilon=1.e-10;

  for(k=0;k<size-1;k++) {
    dx[k]=xp[k+1]-xp[k];
    dy[k]=yp[k+1]-yp[k];
    dn[k]=hypot(dx[k],dy[k]);
    if(dn[k]==0) {
      printf("colocated points in polygon : %d %d at %lf %lf\n",k,k+1, xp[k], yp[k]);
      colocated=true;
      }
    }

/*------------------------------------------------------------------------------
  seek a convenient arbitrary direction*/
  ang=0.0;
  tx=1.0;
  ty=0.0;
  nx=-ty;
  ny=+tx;

  for(l=0;l<size-1;l++) {
    delta[l]=nx*dx[l]+ny*dy[l];
    if (fabs(delta[l]) < error*dn[l]) {
      ang=ang+0.0001;
      if(ang > M_PI+M_PI/4.) {
        STDOUT_BASE_LINE("unsafe polygone, cannot find a scanning direction: %d %d %lf %lf error=%lf\n",l,size,dn[l],delta[l],error);
        ang=0.0;
        error/=2.0;
        }
      tx=cos(ang);
      ty=sin(ang);
      nx=-ty;
      ny=+tx;
      l=-1;
      }
    }
 
/*------------------------------------------------------------------------------
  screen the polygone's segment*/
  for(l=0;l<size-1;l++) {
      ux=lon-xp[l];
      uy=lat-yp[l];
      alpha=(nx*ux+ny*uy)/delta[l];
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      (P,t) is the arbitrary testing half line
      AB is one polygone segment
      B' is the ortogonal (along t) projection on (A,n) axis
      P is (lon,lat)
      P' is the ortogonal (along t) projection on (A,n) axis
      I is the intersection point between (P,t) and (A,AB)

      delta[l] is the B and B'  points' coordinate in (A,n) axis
      nx*ux+ny*uy is the P point coordinate in (A,n) axis
      alpha=(nx*ux+ny*uy)/delta[l] is the barycentric coordinate of the
      intersection point P' in (A,AB') axis
      There is an intersection if 0< alpha < 1
      
xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
/*
      if (alpha == 0.) {
        pos=plg_point_boundary;
        retun(pos);
        }
*/
        if(check && fabs(alpha*(1.-alpha))<1.e-02) {
          printf("check #1: %lf %lf %lf %lf %lf %lf %d %d \n",lon,lat,alpha,beta,dn[l],delta[l],l,*in);
          }

      if ((alpha >= 0.0)&&(alpha < 1.0)) {
/*------------------------------------------------------------------------------
        there is an intersection*/
/*------------------------------------------------------------------------------
        x,y is the intersection position*/
//         x=xp[l]+alpha*dx[l];
//         y=yp[l]+alpha*dy[l];
//         vx=x-lon;
//         vy=y-lat;
        vx=xp[l]-lon;
        vy=yp[l]-lat;
        vx+=alpha*dx[l];
        vy+=alpha*dy[l];
/* *------------------------------------------------------------------------------
        beta is the barycentric coordinate of the intersection point in (A,t) axis */
        beta=vx*tx+vy*ty;
        if(fabs(beta)<epsilon) {
          pos=PLG_POINT_BOUNDARY;
          goto end;
          }
        if(beta>0.) {
          *in=(!(*in));
          }
        if(check) {
          printf("check #1: %lf %lf %lf %lf %lf %lf %d %d \n",lon,lat,alpha,beta,dn[l],delta[l],l,*in);
          }
        }
     }

  if ((*in)==1) {
    pos=PLG_POINT_INTERIOR;
    }
  else {
    pos=PLG_POINT_EXTERIOR;
    }

 end:
  
  delete[] dx;
  delete[] dy;
  delete[] dn;
  delete[] delta;

  return(pos);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_TestInterior(double lon, double lat, plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///test whether a point is on a polygone
/**
\param lon
\param lat
\param *polygones array[npolygones] of plg_t
\param npolygones see above
\return PLG_POINT_BOUNDARY if on a polygone
*/
/*----------------------------------------------------------------------------*/
{
  int tmp,in=0;
  int k;

  /** call plg_single() for each element of \c polygones */
  for(k=0;k<npolygones;k++) {
    tmp=plg_single( polygones[k], lon, lat,&in);
    if(tmp==PLG_POINT_BOUNDARY) break;
    }
  return(tmp);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_TestInterior(double lon, double lat, vector<plg_t> polygons, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///test whether a point is on a polygone
/**
\param lon
\param lat
\param polygons vector of plg_t
\return PLG_POINT_BOUNDARY if on a polygone
*/
/*----------------------------------------------------------------------------*/
{
  int tmp,in=0;
  int k;

  /** call plg_single() for each element of \c polygones */
  for(k=0;k<polygons.size();k++) {
    tmp=plg_single( polygons[k], lon, lat, &in, mode);
    if(tmp==PLG_POINT_BOUNDARY) break;
    }
  return(tmp);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_TestInterior(const plg_t & polygon,const double *lon,const double *lat, int count, char *in, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int size;
  int k,n,lmin;

  double *dx=NULL,*dy=NULL,*dn=NULL,*delta=NULL;
  double epsilon,error,tx,ty,ang,nx,ny;
  double cmin;

  size=polygon.npt;

  exitIfNull(
    dx=new double[size]
    );
  exitIfNull(
    dy=new double[size]
    );
  exitIfNull(
    dn=new double[size]
    );
  exitIfNull(
    delta=new double[size]
    );

  error=1e-01;
  error=1.e-04;

  epsilon=1.e-12;

  ang=M_PI/4.-0.001;
  cmin=1.e+10;
  tx=cos(ang);
  ty=sin(ang);
  nx=-ty;
  ny=+tx;

  for(k=0;k<size-1;k++) {
    dx[k]=polygon.x[k+1]-polygon.x[k];
    dy[k]=polygon.y[k+1]-polygon.y[k];
    dn[k]=sqrt(dx[k]*dx[k]+dy[k]*dy[k]);
    if((dx[k]==0.)&&(dy[k]==0.)) {
      printf("unsafe polygone, 2 successive points colocated (%d %d, C-numbering)\n",k,k+1);
      printf("eliminating the redundant point will make it safe\n");
      return(-1);
      }
    }
  nx=-ty;
  ny=+tx;

  for(int l=0;l<size-1;l++) {
    delta[l]=nx*dx[l]+ny*dy[l];
    if (fabs(delta[l]) < error*dn[l]) {
      if(fabs(delta[l])/dn[l] < cmin ) {
        cmin=fabs(delta[l])/dn[l];
        lmin=l;
        }
      ang=ang+0.001;
      if(ang > M_PI+M_PI/4.) {
        STDOUT_BASE_LINE("unsafe polygone, cannot find a scanning direction: %d %d %lf %lf \n",l,size,dn[l],delta[l]);
        exit(-1);
        }
      tx=cos(ang);
      ty=sin(ang);
      nx=-ty;
      ny=+tx;
      cmin=1.e+10;
      l=-1;
      }
    }

  if(debug)  printf("polygone: size=%d, seeking axis=%lf degrees\n",size,ang*r2d);
  int nprocs=initialize_OPENMP(-1, 0);
#pragma omp parallel for if(nprocs>1)
  for(n=0;n<count;n++) {
    double ux,uy,x,y,vx,vy;
    double alpha,beta;
    int pos;
    if(in[n]==-1) continue;
    for(int l=0;l<size-1;l++) {
      ux=lon[n]-polygon.x[l];
      uy=lat[n]-polygon.y[l];
      alpha=(nx*ux+ny*uy)/delta[l];
      if (alpha == 0.) {
        pos=PLG_POINT_BOUNDARY;
        printf("intersection at point %d, node %d \n",n,l);
        }
      else if (alpha == 1.) {
        pos=PLG_POINT_BOUNDARY;
        printf("intersection at point %d, node %d \n",n,l);
        }
/*------------------------------------------------------------------------------
      there is an intersection*/
      else if ((alpha >= 0.0)&&(alpha < 1.0)) {
/*------------------------------------------------------------------------------
        x,y is the intersection position*/
        x=polygon.x[l]+alpha*dx[l];
        y=polygon.y[l]+alpha*dy[l];
        vx=x-lon[n];
        vy=y-lat[n];
        beta=vx*tx+vy*ty;
        if(fabs(beta)<epsilon) {
          in[n]=-1;
          break;
          }
        if(beta>0.) {
          in[n]=!in[n];
          }
        }
      }
    }

  delete[] dx;
  delete[] dy;
  delete[] dn;
  delete[] delta;

//   for(n=0;n<count;n++) {
//     switch (in[n]) {
//       case -1:
//         in[n]=PLG_POINT_BOUNDARY;
//         break;
//
//       case 0:
//         in[n]=PLG_POINT_EXTERIOR;
//         break;
//
//       case 1:
//         in[n]=PLG_POINT_INTERIOR;
//         break;
//
//       }
//     }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<class T> char *plg_TestInterior_template(const double *lon,const double *lat, int count, const T & polygons, int npolygons, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int k;
  char *position=NULL;
  int granularity=100;
  int milestone, goal;
  
  if(npolygons<100) granularity=min(20, npolygons);
  
  milestone=npolygons/granularity;
  goal=granularity*milestone;

  position=aset(count,'\0');
  
  for(k=0;k<npolygons;k++) {
    if(verbose>0){
      if(k % milestone ==0) printf("%3d percent done (%d/%d)\n", 100*(k/milestone)*milestone/goal,k,npolygons);
      printf("treating polygon %d in %d\n",k,npolygons);
      }
    status=plg_TestInterior( polygons[k], lon, lat, count, position,debug);
    if(status!=0) {
      STDOUT_BASE_LINE("polygon %d found to be unsafe, please correct and retry\n",k);
      exit(-1);
      }
    }

  for(n=0;n<count;n++) {
    switch (position[n]) {
      case -1:
        position[n]=PLG_POINT_BOUNDARY;
        break;

      case 0:
        position[n]=PLG_POINT_EXTERIOR;
        break;

      case 1:
        position[n]=PLG_POINT_INTERIOR;
        break;

      }
    }

  return(position);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *plg_TestInterior(const double *lon,const double *lat, int count,const plg_t *polygones, int npolygones, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *position;
  
  position=plg_TestInterior_template(lon,lat,count,polygones,npolygones,verbose, debug);
  
  return position;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *plg_TestInterior(const double *lon,const double *lat, int count,const vector<plg_t> & polygons, int verbose, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *position;
  
  position=plg_TestInterior_template(lon,lat,count,polygons,polygons.size(),verbose, debug);
  
  return position;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_GridInterior_atom(const plg_t & polygon,const grid_t & grid, char *in)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  optimisation may not work on curvilinear grids

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  int size;
  int k,l,lmin;
  double *dx=NULL,*dy=NULL,*dn=NULL,*delta=NULL;
  double epsilon,error,tx,ty,ang,nx,ny;
  double cmin;
  frame_t frame;

  size=polygon.npt;
  
  if(grid.modeH!=2) TRAP_ERR_EXIT(ENOEXEC,"Programming error: %s(,grid,) must be called with grid.modeH==2\n",__func__);

  exitIfNull(
    dx=new double[size]
    );
  exitIfNull(
    dy=new double[size]
    );
  exitIfNull(
    dn=new double[size]
    );
  exitIfNull(
    delta=new double[size]
    );

//  error=1e-03;
  error=1.e-05;

  epsilon=1.e-12;

  ang=M_PI/4.-0.001;
  cmin=1.e+10;
  tx=cos(ang);
  ty=sin(ang);
  nx=-ty;
  ny=+tx;

  for(k=0;k<size-1;k++) {
    dx[k]=polygon.x[k+1]-polygon.x[k];
    dy[k]=polygon.y[k+1]-polygon.y[k];
    dn[k]=sqrt(dx[k]*dx[k]+dy[k]*dy[k]);
    if((dx[k]==0.)&&(dy[k]==0.)) {
      printf("duplicated consecutive points in polygon: %d %d\n",k,k+1);
      printf("position: lon=%lf lat=%lf x=%lf y=%lf\n",polygon.t[k],polygon.p[k],polygon.x[k],polygon.y[k]);
      printf("position: lon=%lf lat=%lf x=%lf y=%lf\n",polygon.t[k+1],polygon.p[k+1],polygon.x[k+1],polygon.y[k+1]);
      }
    }
  nx=-ty;
  ny=+tx;

  for(l=0;l<size-1;l++) {
    delta[l]=nx*dx[l]+ny*dy[l];
    if (fabs(delta[l]) < error*dn[l]) {
      if(fabs(delta[l])/dn[l] < cmin ) {
        cmin=fabs(delta[l])/dn[l];
        lmin=l;
        }
      ang=ang+0.001;
      if(ang > M_PI+M_PI/4.) {
        STDOUT_BASE_LINE("erreur de polygone: %d %d %lf %lf \n",l,size,dn[l],delta[l]);
        exit(-1);
        }
      tx=cos(ang);
      ty=sin(ang);
      nx=-ty;
      ny=+tx;
      cmin=1.e+10;
      l=-1;
      }
    }

  frame=plg_cartesian_minmax(&polygon, 1);

//  printf("polygone: size=%d, seeking axis=%lf degrees\n",size,ang*r2d);

  for(l=0;l<size-1;l++) {
    const double dd=delta[l];
    const double xx=polygon.x[l];
    const double yy=polygon.y[l];
#pragma omp parallel for
    for(int j=0;j<grid.ny;j++) {
//     if( (lon[n]<frame.xmin) || (lon[n]>frame.xmax)  || (lat[n]<frame.ymin) || (lat[n]>frame.ymax) ) {
//       continue;
//       }
/*xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx

      Intersection of testing axis with the j-th horizontal grid axis

      x=x0+r tx
      y=y0+r ty

      x=x0+(y-y0) tx/ty

xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx*/
      double y=grid.y[j*grid.nx];
      double ratio=tx/ty;
      double x1=polygon.x[l]  +(y-polygon.y[l]  )*ratio;
      double x2=polygon.x[l+1]+(y-polygon.y[l+1])*ratio;
      int i1=(int) min(floor((x1-grid.xmin)/grid.dx), floor((x2-grid.xmin)/grid.dx))-2;
      int i2=(int) max( ceil((x1-grid.xmin)/grid.dx),  ceil((x2-grid.xmin)/grid.dx))+2;
      i1=max(0,i1);
      i2=min(i2,grid.nx);
      for(int i=i1;i<i2;i++) {
        double x,ux,uy,vx,vy;
        int pos;
        size_t n=(size_t) j*(size_t) grid.nx+(size_t) i;
//         ux=grid.x[n]-polygon.x[l];
//         uy=grid.y[n]-polygon.y[l];
        ux=grid.x[n]-xx;
        uy=grid.y[n]-yy;
/*------------------------------------------------------------------------------
        alpha is barycentric coordinate of intersection point along segment  */
//         alpha=(nx*ux+ny*uy)/delta[l];
        const double alpha=(nx*ux+ny*uy)/dd;
        if (alpha == 0.) {
          pos=PLG_POINT_BOUNDARY;
          printf("intersection au noeud: %d %d \n",n,l);
          }
        else if (alpha == 1.) {
          pos=PLG_POINT_BOUNDARY;
          printf("intersection au noeud: %d %d \n",n,l);
          }
/*------------------------------------------------------------------------------
        there is an intersection*/
        else if ((alpha > 0.0) and (alpha < 1.0)) {
/*------------------------------------------------------------------------------
          x,y is the intersection position*/
          x=polygon.x[l]+alpha*dx[l];
          y=polygon.y[l]+alpha*dy[l];
          vx=x-grid.x[n];
          vy=y-grid.y[n];
/*------------------------------------------------------------------------------
          beta is barycetric coordinate of intersection point across segment*/
          const double beta=vx*tx+vy*ty;
          if(fabs(beta)<epsilon) {
            in[n]=-1;
            break;
            }
          if(beta>0.) {
            in[n]=(!(in[n]));
            }
          }
        }
      }
    }

  delete[] dx;
  delete[] dy;
  delete[] dn;
  delete[] delta;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<class T> char *plg_GridInterior_template(const grid_t & grid, const T & polygons, int npolygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  int k;
  char *position=NULL;
  int count;

  int nRequestedProcs=-1;
  int nprocs __attribute__((unused)) =initialize_OPENMP(nRequestedProcs, 0);
  
  count=grid.nx*grid.ny;

  position=aset(count,(char)0);
  
  for(k=0;k<npolygons;k++) {
//     printf("treating polygone %d\n",k);
    status=plg_GridInterior_atom( polygons[k], grid, position);  //  < need optimisation
    }

  for(n=0;n<count;n++) {
    switch (position[n]) {
      case -1:
        position[n]=PLG_POINT_BOUNDARY;
        break;

      case 0:
        position[n]=PLG_POINT_EXTERIOR;
        break;

      case 1:
        position[n]=PLG_POINT_INTERIOR;
        break;

      }
    }

  return(position);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *plg_test_grid(const grid_t & grid,const plg_t *polygons, int npolygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *position;
  
  position=plg_GridInterior_template(grid,polygons,npolygons);
  
  return position;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  char *plg_test_grid(const grid_t & grid, const vector<plg_t> & polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  char *position;
  
  position=plg_GridInterior_template(grid,polygons,polygons.size());
  
  return position;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_secantpoint(double x1[2],double y1[2],double x2[2],double y2[2],double *x,double *y,double *a,double *b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  first  line equation : M=P+a U
  second line equation : N=Q+b V

  intersection for M=N : a U - b V = Q-P = W
  
  a Ux -b Vx = (Q-P)x
  a Uy -b Vy = (Q-P)x
  
  Ux -Vx
  Uy -Vy      determinant = Uy Vx - Ux Vy = - U * V
  
  a = (W * V)/(U * V)
  b = (U * W)/(U * V)

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  double epsilon=0.;
  double u[2],v[2],a_,b_,dx,dy,d,xx,yy;
  
  if(a==0) a=&a_;
  if(b==0) b=&b_;

  u[0]=x1[1]-x1[0];
  u[1]=y1[1]-y1[0];

  v[0]=x2[1]-x2[0];
  v[1]=y2[1]-y2[0];

  dx=x2[0]-x1[0]; // Wx
  dy=y2[0]-y1[0]; // Wy

  d=u[1]*v[0]-u[0]*v[1];

  if(!isnormal(d)) {
/* *----------------------------------------------------------------------------
    lines are parallel or ill-defined */
    if(dx==0 && dy==0) return(PLG_LINES_NOT_SECANT);
    else return(PLG_LINES_NOT_SECANT);
    }

  *a=( v[0]*dy-v[1]*dx)/d;
  *b=( u[0]*dy-u[1]*dx)/d;

  if(*a < 0.) {
    return(PLG_LINES_NOT_SECANT);
    }
  if(*a > 1.) {
    return(PLG_LINES_NOT_SECANT);
    }
  if(*b < 0.) {
    return(PLG_LINES_NOT_SECANT);
    }
  if(*b > 1.) {
    return(PLG_LINES_NOT_SECANT);
    }

  xx=x1[0]+*a*u[0];
  yy=y1[0]+*a*u[1];
  
  dx=x1[0]-xx;
  dy=y1[0]-yy;
  if(fabs(dx)+fabs(dy)==0) *a=0.0;
  
  dx=x1[1]-xx;
  dy=y1[1]-yy;
  if(fabs(dx)+fabs(dy)==0) *a=1.;
  
  dx=x2[0]-xx;
  dy=y2[0]-yy;
  if(fabs(dx)+fabs(dy)==0) *b=0.0;
  
  dx=x2[1]-xx;
  dy=y2[1]-yy;
  if(fabs(dx)+fabs(dy)==0) *b=1.;
  
  *x=xx;
  *y=yy;
  
  dx=u[1]*x1[0]-u[0]*y1[0];
  dy=v[1]*x2[0]-v[0]*y2[0];
  
  xx=-(+u[0]*dy-v[0]*dx)/d;
  yy=-( u[1]*dy-v[1]*dx)/d;

  if(*a == 0.) {
    return(PLG_LINES_SECANT_AT_EXTRIMITY);
    }
  if(*a == 1.) {
    return(PLG_LINES_SECANT_AT_EXTRIMITY);
    }

  if(*b <= epsilon) {
    return(PLG_LINES_SECANT_AT_EXTRIMITY);
    }
  if(*b >= 1.-epsilon) {
    return(PLG_LINES_SECANT_AT_EXTRIMITY);
    }
  
  return(PLG_LINES_SECANT);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_secantpoint(const line_t & A, const line_t & B, double *x,double *y, double *a, double *b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int flag;
  
  double x1[2], y1[2], x2[2], y2[2];
  double a_,b_;
  if(a==0) a=&a_;
  if(b==0) b=&b_;

  
  x1[0]=A.point[0].t;
  x1[1]=A.point[1].t;
  x2[0]=B.point[0].t;
  x2[1]=B.point[1].t;
  
  y1[0]=A.point[0].p;
  y1[1]=A.point[1].p;
  y2[0]=B.point[0].p;
  y2[1]=B.point[1].p;
  
  flag=plg_secantpoint(x1, y1, x2, y2, x, y, a, b);
  
  return flag;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_secantpoint(const line_t & A, const line_t & B, plg_point_t & point, double *a, double *b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int flag;
  
  double x,y;
    
  flag=plg_secantpoint(A, B, &x, &y, a, b);
  
  point.t=x;
  point.p=y;
 
  point.x=x;
  point.y=y;
  
  return flag;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_secantpoint(const line_t & A, const line_t & B, point2D_t & point, double *a, double *b)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int flag;
  double x,y;
    
  flag=plg_secantpoint(A, B, &x, &y, a, b);
  
//   point.t=x;
//   point.p=y;
 
  point.x=x;
  point.y=y;
  
  return flag;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_secantpoints(const plg_t & p, const plg_t & q, point2D_t **points, double **angles, double **aindexes, double **bindexes)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  scan intersections between 2 polygons

 **points pointer to array    : if not NULL, as is the default, will be allocated
                                to contain the coordinates
 **angles pointer to array    : if not NULL, as is the default, will be allocated
                                to contain the angles
 **aindexes pointer to array  : if not NULL, as is the default, will be allocated
                                to contain the indexes of the intersected segments in p
 **bindexes pointer to array  : if not NULL, as is the default, will be allocated
                                to contain the indexes of the intersected segments in q

  return the number of intersected points
  
  based on (damien-style) 2 passes loop (first count and allocate, second registers data

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  point2D_t *points0=0;
  int countMade,npoints,k,l,status;
  double xa[2], ya[2], xb[2], yb[2], x, y, ka, kb;
  
  for(countMade=0;;countMade++){
    
    npoints=-1;
    
    for(k=0;k<p.npt-1;k++){
      memcpy(xa,&p.x[k],2*sizeof(double));
      memcpy(ya,&p.y[k],2*sizeof(double));
      for(l=0;l<q.npt-1;l++){
        memcpy(xb,&q.x[l],2*sizeof(double));
        memcpy(yb,&q.y[l],2*sizeof(double));
        
        status=plg_secantpoint(xa, ya, xb, yb, &x, &y, &ka, &kb);
        
        if(status==PLG_LINES_SECANT_AT_EXTRIMITY && kb<=0.)
          status=PLG_LINES_SECANT;
        
        if(status!=PLG_LINES_SECANT) continue;
        
        npoints++;
/*------------------------------------------------------------------------------
        if first pass, just make count */
        if(countMade==0) continue;
        
/*------------------------------------------------------------------------------
        second pass, register information */
        points0[npoints].init(x,y);
        if(angles!=0){
          const vector2D_t
            a(xa[0]-xa[1],ya[0]-ya[1]),
            b(xb[0]-xb[1],yb[0]-yb[1]);
          (*angles)[npoints]=a^b;
          }
        if(aindexes!=0)
          (*aindexes)[npoints]=k+ka;
        if(bindexes!=0)
          (*bindexes)[npoints]=l+kb;
        }
      }
    
/*------------------------------------------------------------------------------
    if second pass, leave */
    if(countMade==1) break;
/*------------------------------------------------------------------------------
    if no information to be returned, leave */
    if(points==0 and angles==0 and bindexes==0) break;
    
/*------------------------------------------------------------------------------
    if no intersections, leave */
    if(npoints<0)return 0;
    
/*------------------------------------------------------------------------------
    end of first pass, allocate arrays if convenient */
    npoints++;
    points0=new point2D_t[npoints];
    if(angles!=0)
      *angles=new double[npoints];
    if(aindexes!=0)
      *aindexes=new double[npoints];
    if(bindexes!=0)
      *bindexes=new double[npoints];
    }
  
  if(points!=0)
    *points=points0;
  
  return (npoints+1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_secantpoints(const vector<plg_t> & p, const line_t line, vector<point2D_t> & points)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  point2D_t point;
  int npoints,k,status;
  double a, b;
  vector<point2D_t> tmp;
  vector<double> coordinate;
  size_t *pos=0;
  
  for(size_t s=0;s<p.size();s++){
    for(k=0;k<p[s].npt-1;k++){
      line_t test(p[s],k,PLG_SPHERICAL);
      status=plg_secantpoint(line, test, point,&a,&b);
      if(status==PLG_LINES_SECANT_AT_EXTRIMITY && b<=0.)
        status=PLG_LINES_SECANT;
      if(status!=PLG_LINES_SECANT) continue;
      coordinate.push_back(a);
      tmp.push_back(point);
      }
    }
  
  pos=sort(coordinate);
  for(size_t k=0; k< tmp.size();k++) {
    points.push_back(tmp[pos[k]]);
    }

  tmp.clear();
  coordinate.clear();
  
  if(pos!=0) delete[]pos;
  
  npoints=points.size();
  return (npoints);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_checkAutoSecant(plg_t *polygones, int npolygones, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double x,y;
  
  vector<int> *flags;
  
  flags=new vector<int>[npolygones];

/*------------------------------------------------------------------------------
  check intersection in a given polygon */
  for(size_t p=0;p<npolygones;p++) {
    
    if(polygones[p].npt==0) continue;
    
    for(size_t k=0;k<polygones[p].npt-1;k++) {
      line_t a=line_t(polygones[p],k,mode);
      for(size_t l=k+2;l<polygones[p].npt-1;l++) {
        if(fabs(k-l)<2) continue;
        line_t b=line_t(polygones[p],l,mode);
        int flag=plg_secantpoint(a,b,&x,&y);
        if(flag==PLG_LINES_SECANT) {
          printf("polygon/segment %d is auto-secant, lines are %d %d\n",p,k,l);
          printf("intersection position: %lfE %lfN\n",x,y);
          flags[p].push_back(k);
          break;
          }
        }
      }
    }

  status=0;
  for(size_t p=0;p<npolygones;p++) {
    if(flags[p].size()!=0) {
      status=-1;
      }
    }
  
  delete [] flags;
  
  if (status!=0) return(-1);
  
  flags=new vector<int>[npolygones];

/*------------------------------------------------------------------------------
  check duplicate points polygon */
#pragma omp parallel for
  for(size_t p=0;p<npolygones;p++) {
    
    if(polygones[p].npt==0) continue;
    
    for(size_t k=0;k<polygones[p].npt-1;k++) {
      for(size_t l=k+1;l<polygones[p].npt-1;l++) {
        double d=plg_distance(polygones[p], k, l, PLG_CARTESIAN);
        int flag=(d==0.0);
        if(flag==1) {
          printf("polygon %d has duplicated points, %d %d\n",p,k,l);
          printf("duplicated position: %lfE %lfN\n",polygones[p].t[k],polygones[p].p[k]);
          flags[p].push_back(k);
          break;
          }
        }
      }
    }

  status=0;
  for(size_t p=0;p<npolygones;p++) {
    if(flags[p].size()!=0) {
      status=-1;
      }
    }
  
  delete [] flags;
  
  if (status!=0) return(-1);
  
  flags=new vector<int>[npolygones];

/*------------------------------------------------------------------------------
  check intersection with a another polygon */
#pragma omp parallel for
  for(size_t p=0;p<npolygones;p++) {
    
    if(polygones[p].npt==0) continue;
    
    for(size_t k=0;k<polygones[p].npt-1;k++) {
      line_t a=line_t(polygones[p],k,mode);
      for(size_t q=p+1;q<npolygones;q++) {
    
        if(polygones[q].npt==0) continue;
    
        for(size_t l=0;l<polygones[q].npt-1;l++) {
          line_t b=line_t(polygones[q],l,mode);
          int flag=plg_secantpoint(a,b,&x,&y);
          if(flag==PLG_LINES_SECANT) {
            printf("polygon/segment %d is secant with %d, lines are %d %d\n",p,q,k,l);
            printf("intersection position: %lfE %lfN\n",x,y);
            flags[p].push_back(q);
            break;
            }
          }
        }
      }
    }
  
  status=0;
  for(size_t p=0;p<npolygones;p++) {
    if(flags[p].size()!=0) {
      status=-1;
      }
    }

  delete [] flags;
    
  return(status);
  
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_checkAutoSecant(vector<plg_t> polygons, int mode)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  plg_array_t tmp;

  tmp=plg_vector2array(polygons);
  status=plg_checkAutoSecant(tmp.p,tmp.n, mode);
  
  delete[] tmp.p;
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int plg_checkSecant(const vector<plg_t> & p, const vector<plg_t> & q, bool debug)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int npoints, global=0;
  point2D_t **points=0;
  double **angles=0;
  double *aindexes=0, *bindexes=0;

/* *----------------------------------------------------------------------------
  check intersection with a another polygon */
// #pragma omp parallel for private(status,npoints)
  for(size_t l=0;l<p.size();l++) {
    
    if(p[l].npt==0) continue;
    
    for(size_t k=0;k<q.size();k++) {
      if(q[k].npt==0) continue;
      npoints=plg_secantpoints(p[l], q[k],  points, angles, &aindexes, &bindexes);
      if(npoints!=0) {
        global+=npoints;
        if(debug) {
          printf("polygons %d %d have %d points of intersection\n",l,k,npoints);
          for(int n=0;n<npoints;n++) {
            int ll=floor(aindexes[n]+0.5);
//             int kk=floor(bindexes[n]+0.5);
            printf("%d %lf %lf, t=%lf p=%lf \n",n,aindexes[n],bindexes[n], p[l].t[ll], p[l].p[ll]);
            }
          }
        }
      }
    }
  
  return(global);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int enforce_landmask_template(const grid_t & grid, plg_t *polygones, int npolygones, T* mask, T value, bool tuned)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  
  set landmask to max(value, mask) at grid closest nodes to a shoreline node
  
  it is aimed to ensure laguna and peninsula proper treatments;
  !!! quite rough, needs to be revised
  
  its major drawback is a shifting effect when the working grid granularity
  is close to the sampling resolution
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
  if(npolygones ==0) return(0);
  
  for(size_t p=0;p<npolygones;p++) {
    
    if(polygones[p].npt==0) continue;
     
    double L=plg_length(polygones[p], PLG_CARTESIAN);
    
    if(L<100*grid.dx) continue;

    if(tuned) {
      plg_t q=plg_resample(polygones[p],grid.dx/3.0, false);
    
      if(q.npt==0) continue;
    
      for(size_t k=0;k<q.npt;k++) {
        const double x=q.x[k];
        const double y=q.y[k];
        int    i,j,m;
        i=NINT((x-grid.xmin)/grid.dx);
        j=NINT((y-grid.ymin)/grid.dy);
        m=grid.nx*j+i;
        mask[m]=max(value,mask[m]);
        }
      q.destroy();
      }
    else {
      for(size_t k=0;k<polygones[p].npt;k++) {
        const double x=polygones[p].x[k];
        const double y=polygones[p].y[k];
        int    i,j,m;
        i=NINT((x-grid.xmin)/grid.dx);
        j=NINT((y-grid.ymin)/grid.dy);
        m=grid.nx*j+i;
        mask[m]=max(value,mask[m]);
        }
      }
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

 int enforce_landmask(const grid_t & grid, plg_t *polygones, int npolygones, float* mask, float value, bool tuned)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=enforce_landmask_template(grid, polygones, npolygones, mask, value, tuned);
  
  return(status);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_landmask01(const grid_t & grid, plg_t *polygones, int npolygones)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  *mask;
  char *interior;
  const int nxy=grid.nx*grid.ny;

  if(npolygones ==0) return(0);

  mask=new float[grid.nx*grid.ny];

//   int nRequestedProcs=-1;
//   int nprocs=initialize_OPENMP(nRequestedProcs);
//
// #pragma omp parallel for private(status) if(nprocs>1)
//   for(j=0;j<grid.ny;j++) {
//     for(int i=0;i<grid.nx;i++) {
//       size_t k=grid.nx*j+i;
//       double lon=grid.x[k];
//       double lat=grid.y[k];
//       mask[k]=plg_TestInterior(lon,lat,polygones,npolygones);
//       }
//     }
    
  interior=plg_test_grid(grid,polygones,npolygones);
  
//   for(size_t p=0;p<npolygones;p++) {
//     
//     if(polygones[p].npt==0) continue;
//     
//     for(size_t k=0;k<polygones[p].npt;k++) {
//       const double x=polygones[p].x[k];
//       const double y=polygones[p].y[k];
//       int    i,j,m;
//       i=NINT((x-grid.xmin)/grid.dx);
//       j=NINT((y-grid.ymin)/grid.dy);
//       m=grid.nx*j+i;
//       interior[m]=1;
//       }
//     }

/*------------------------------------------------------------------------------
  fixes mask and interior incompatible types */  
  valcpy(mask,interior,nxy);
  
  delete[] interior;

  return(mask);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_landmask01(const grid_t & grid, vector<plg_t> polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  *mask;
  char *interior;
  const int nxy=grid.nx*grid.ny;

  if(polygons.size() ==0) return(0);

  mask=new float[grid.nx*grid.ny];
    
  interior=plg_test_grid(grid,polygons);
  
/*------------------------------------------------------------------------------
  fixes mask and interior incompatible types */  
  valcpy(mask,interior,nxy);
  
  delete[] interior;

  return(mask);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *set_landmask02(const grid_t & grid, vector<plg_t> & polygons, bool checks)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    status;
  float  *mask=NULL;
  char   *interior;
  int weird=0;
  frame_t frame=frame_t(grid.xmin,grid.xmax,grid.ymin,grid.ymax);
  const int nxy=grid.nx*grid.ny;
  bool repair=true, verbose=false, debug=false;;
  
  if(polygons.size() ==0)
    return 0;
  
  if(checks) {
    weird=0;
    for(int k=0;k<polygons.size();k++) {
      if(polygons[k].npt<3) {
        printf("checking %4dth polygon, id=%d: polygon has only 2 or less points\n",k, polygons[k].id);
        printf("first position: %lfE %lfN\n",polygons[k].t[0],polygons[k].p[0]);
        printf("%6d %d\n", polygons[k].id, polygons[k].npt);
        weird++;
        }
      }
    for(int k=0;k<polygons.size();k++) {
      if(plg_isclosed(polygons[k])==0) {
        printf("checking %4dth polygon, id=%d: polygon not closed\n",k, polygons[k].id);
        printf("first position: %lfE %lfN\n",polygons[k].t[0],polygons[k].p[0]);
        printf("%6d %d\n", polygons[k].id, polygons[k].npt);
        weird++;
        }
      }
    if(weird!=0) {
      TRAP_ERR_EXIT(-1, "selection polygons not safe, abort\n");
      }
    status=plg_CheckDuplicated(polygons, repair, (string) "./", verbose);
    if(status!=0) {
      TRAP_ERR_EXIT(-1, "selection polygons not safe, abort\n");
      }
    }
  
  status=plg_recale(polygons, frame, PLG_SPHERICAL);

  interior=plg_TestInterior(grid.x, grid.y,nxy, polygons, 0, debug);

  poc_copy(mask,interior,nxy);

  delete[] interior;

  return(mask);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  float *set_landmask02(const grid_t & grid, plg_t *polygons, int npolygons, bool checks)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  vector<plg_t> p;
  float  *mask=NULL;
  
  TRAP_ERR_EXIT(ENOEXEC,"bug : p is not initialised with polygons\n");
  mask=set_landmask02(grid, p, checks);
  
  p.clear();
  
  return(mask);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

float *set_landmask03(const grid_t & grid, vector<plg_t> polygons)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int    status;
  float  *mask;
  char *interior;

  if(polygons.size() ==0) return(0);

  mask=aset(grid.Hsize(),-1.f);
  status=plg_SortBySize(polygons,0);
  
  for(int p=0;p<polygons.size();p++) {
    interior=plg_test_grid(grid,&polygons[p],1);
    for(int j=0;j<grid.ny;j++) {
      for(int i=0;i<grid.nx;i++) {
        size_t k=grid.nx*j+i;
        if(interior[k]==PLG_POINT_INTERIOR) mask[k]=p;
        }
      }
    delete[] interior;
    }

  return(mask);
}
