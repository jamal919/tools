
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief map interpolation definitions
*/
/*----------------------------------------------------------------------------*/

#include <complex>

#include "map-classes.h" //for grid_t
#include "maths.h" //for vector2D_t
#include "poc-list.hpp" /* for poc_deque_t */

#if TUGO
#include "fe.h" //for findInList
#else
#include "fe-proto.h" //for findInList
#endif


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  bool map_check_circular(const grid_t & grid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{
  bool circular;
  double x0,xn;
  
  if(map_recale(grid,grid.xmax)-grid.xmin<180.)
    return false;
  
  int m0, mn;
  
  mn=grid.nx-1;
  
  switch(grid.modeH){
  
  case 0:
/*------------------------------------------------------------------------------
    patch sordide a revoir */
    x0=grid.xmin;
    xn=grid.xmax+grid.dx;
    xn=degree_recale(xn,x0);
    circular=is_between(xn-grid.dx*.5,x0,xn+grid.dx*.5);
    return circular;
    break;
  
  case 1:
    m0=0;
    break;
  
  case 2:
    m0=(grid.ny/2)*grid.nx;
    break;
  
  case -2:
    m0=grid.nx/2;
    mn=(grid.ny-1)*grid.nx;
    break;
  
  case MODEH_REDUCED_GG:
    return true;
  
  default:
    TRAP_ERR_EXIT(ENOEXEC,"%s programming error: grid_t::modeH=%d\n",__func__,grid.modeH);
    }
  
  mn+=m0;
  
  x0=grid.x[m0];
  xn=grid.x[mn]+grid.dx;
  xn=degree_recale(xn,x0);
  
  circular=is_between(xn-grid.dx*.5,x0,xn+grid.dx*.5);
  
  return circular;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void set_grid_list(grid_t *grid,int wrap,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///set grid_t::list
/**
\param *grid
\param wrap
\param verbose If >=8: testing mode!
*/
/*----------------------------------------------------------------------------*/
{
  int j;
  
  if(grid->nx==1 or grid->ny==1)
    TRAP_ERR_RETURN(,verbose,"%s: fail safe: %dx%d\n",__func__,grid->nx,grid->ny);
  if(grid->xmin>=grid->xmax || grid->ymin>=grid->ymax)
    TRAP_ERR_EXIT(EIO,"%s: data error: [%g;%g][%g;%g]\n",__func__,grid->xmin,grid->xmax,grid->ymin,grid->ymax);
  
/*--------------------------------------------------------------------*//**<h1>
  also wraps grid around </h1>*/
  grid->circular=wrap && map_check_circular(*grid);

  if(grid->circular){
    STDERR_BASE_LINE("*** PATCH FOR WRAPPED AROUND GRID ***\n");
    
    double *x=grid->x,*y=grid->y;
    const int
      nx0=grid->nx,
      ny0=grid->ny;
    int m;
    if(grid->modeH==-2)
      grid->ny++;
    else
      grid->nx++;
    
    grid->xmax=grid->xmin+360.;
    
    switch(grid->modeH){
      case 1:
        grid->x=new double[grid->nx];
        
        valcpy(grid->x,x,nx0);
        grid->x[nx0]=x[0]+360.;
        
        deletep(&x);
        break;
    
      case 2:
      case -2:{
        int nxy=grid->Hsize();
        grid->x=new double[nxy];
        grid->y=new double[nxy];
        
        if(grid->modeH==2){
          int m0;
          for(m=0,m0=0;m<nxy;m+=grid->nx,m0+=nx0){
            valcpy(&grid->x[m],&x[m0],nx0);
            grid->x[m+nx0]=x[m0];
            valcpy(&grid->y[m],&y[m0],nx0);
            grid->y[m+nx0]=y[m0];
            }
          }
        else{
          const int nxy0=nx0*ny0;
          valcpy(grid->x,x,nxy0);
          valcpy(grid->y,y,nxy0);
          valcpy(&grid->x[nxy0],x,nx0);
          valcpy(&grid->y[nxy0],y,nx0);
          }
        
        deletep(&x);
        deletep(&y);
        }break;
    
      case MODEH_REDUCED_GG:
        break;
    
      default:
        TRAP_ERR_EXIT(ENOEXEC,"%s programming error: grid_t::modeH=%d\n",__func__,grid->modeH);
        }
    }
  
/*-----------------------------------------------------------------------------
  finally make list */
  list_t *list;
  grid_t gridList;
  
  if(grid->nx<=0 || grid->ny<=0 || abs(grid->modeH)!=2){
    list=NULL;
    goto checkout;
    }
  
  gridList=*grid;
  gridList.x=NULL;
  gridList.y=NULL;
  gridList.modeH=0;
  gridList.nx=max(grid->nx,2);
  gridList.ny=max(grid->ny,2);
  map_set_dxdy(&gridList);
  if(verbose>1) gridList.print(cerr,__FILE__ ":  ");
  
  list=new list_t[1];
  
  list->init(gridList.nx,gridList.ny);
  list->dx=gridList.dx;
  list->dy=gridList.dy;
  
  bool verboseWrap;
  
#pragma omp parallel for
  for(j=0;j<grid->ny;j++){
    int i,in,jn,mn;
    int m=j*grid->nx;
    for(i=0;i<grid->nx;i++,m++){
      range_t<double> xR,yR;
      
      for(jn=max(j-1,0);jn<=j+1;jn++){
        if(jn>=grid->ny) break;
        mn=jn*grid->nx;
        
        for(in=max(i-1,0);in<=i+1;in++){
          if(in>=grid->nx) break;
          
          if(isinf(xR.min))
            xR << grid->x[mn+in];
          else
            xR << degree_recale(grid->x[mn+in],xR.min);
          yR << grid->y[mn+in];
          }
        }
      
      if(!isfinite(xR) || !isfinite(yR)) continue;
      const double xm=degree_recale(grid->x[m],xR.min);
      const double ym=grid->y[m];
      if(isnan(xm) || isnan(ym)) continue;
      
      xR.min=(xR.min+xm)*.5;
      xR.max=(xR.max+xm)*.5;
      yR.min=(yR.min+ym)*.5;
      yR.max=(yR.max+ym)*.5;
      
      verboseWrap=
        verbose>=6 ||
        (yR.in(-89.,89.) && xR.d()>120.) ||
        (verbose>=2 && j==grid->ny/2 && i==grid->nx/2);
      if(verboseWrap)
        STDERR_BASE_LINE("[%d,%d](%g,%g):fe_UpdateList(,,[%g,%g],[%g,%g],%d)\n",i,j,grid->x[m],grid->y[m],yR.min,yR.max,xR.min,xR.max,m);
      
      fe_UpdateList(gridList,list,yR,xR,m,verboseWrap);
      
      if(verbose>8 && i && j) TRAP_ERR_EXIT(ENOEXEC,"testing\n");
      }
    }
  if(verbose>=4){
    size_t n;
    int i;
    double sum=0.;
    for(j=0;j<list->ny;j++){
      for(i=0;i<list->nx;i++){
        vector<int> &eisij=list->elements[i][j];
        n=eisij.size();
        if(verbose>=4){
          STDERR_BASE_LINE("[%d,%d].size()=%u",i,j,n);
          for(int i=0;i<n;i++) fprintf(stderr," %d",eisij[i]);
          fprintf(stderr,"\n");
          }
        sum+=n;
        }
      }
    STDERR_BASE_LINE_FUNC("mean([,].size())=%g\n",sum/gridList.Hsize());
    }
  if(verbose>=8) TRAP_ERR_EXIT(ENOEXEC,"testing\n");
  
checkout:
  if(grid->list) delete grid->list;
  
  grid->list=list;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  size_t m4buf(const grid_t &grid, size_t m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
  correct index of a 2D grid in the case it is connex
  
  08/10/2018 : connex is not the correct flag, circular is
  
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc*/
{
//   if(not grid.connex)
//     return m;
  if(not grid.circular)
    return m;
  
  int i,j;
  grid.ij(m,&i,&j);
  
  if(grid.modeH==-2){
    if(j==grid.ny-1)
      return i;
    }
  else{
    const size_t nx0=grid.nx-1;
    
    if(i==nx0)
      i=0;
    
    m=nx0*j+i;
    }
  
  return m;
}


template<typename T> T bilinear_interpolation(const grid_t &grid,double lon,double lat,int m,const T *buf,const T & mask,int verbose);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int index_interpolation02(const grid_t &grid,double t,double p,const list_t & list,const T *buf,const T & mask,T *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///gives indexes of a 2D grid
/**
\return index of closest
*/
/*----------------------------------------------------------------------------*/
{
  if(abs(grid.modeH)!=2)TRAP_ERR_EXIT(ENOEXEC,"%s should not have been called with grid.modeH=%d\n",__func__,grid.modeH);
  
  double x,y;
  const vector<int> *found;
  
  found=findInList(list,grid,list.dx,list.dy,t,p,&x,&y,verbose);
  if(found==0){
    if(z!=0)
      *z=mask;
    return -1;
    }
  
  poc_deque_t<int> eis(*found);
  int i,imin;
  int64_t m,m_,mmin;
  double d,dmin;
  T zm,z0=mask;
  
  while(z0==mask){
    if(verbose)
      STDERR_BASE_LINE_FUNC("(%g,%g)%d\n",t,p,eis.size());
    
    mmin=-1;
    dmin=+INFINITY;
    
    for(i=0;i<eis.size();i++){
      m=eis[i];
      m_=m4buf(grid,m);
      zm=buf[m_];
      if(zm==mask) continue;
      
      d=hypot(x-degree_recale(grid.x[m],x),y-grid.y[m]);
      if(verbose>0 or isinf(d)) {
        STDERR_BASE_LINE("(%g,%g)%d:%g,",grid.x[m],grid.y[m],m,d);
        cerr << zm << "\n";
        }
      if(dmin>d){
        dmin=d;
        mmin=m;
        imin=i;
        }
      }

    if(mmin<0){
      if(z!=0)
        *z=mask;
      return -2;
      }
    
    z0=bilinear_interpolation(grid,t,p,mmin,buf,mask,verbose);
    if(verbose) {
      STDERR_BASE_LINE_FUNC("(%g,%g)",t,p);
      cerr << z0 << "\n";
      }
    
    if(z0==mask) {
      eis.erase(imin);
      }
    }
  
  if(z!=0)
    *z=z0;
  return mmin;
}


template<typename T> void incwz(T zm,T mask,double wm,double *wt,double *w,T *z);


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T bilinear_reduced_gg_interpolation(const grid_t &grid,double lon,double lat,int xi,int i,int j,const T *buf,const T & mask,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///bilinear part of interpolation
/**
\param lon longitude in degrees
\param lat latitude in degrees
*/
/*----------------------------------------------------------------------------*/
{
  int status;
  int j0,j1,xi0,xi1,i0,i1;
  
  if(verbose>0)
    STDERR_BASE_LINE("\n");
  
  i0=i1=i;
  
  if(lat>grid.y[j] or j==grid.ny-1){/* NOTE: latitude are decreasing */
    j0=j-1;
    j1=j;
    xi0=grid.reduced_xi[j0];
    status=ind1D(&grid.x[xi0],grid.reduced_nx[j0],lon,&i0);
    xi1=xi;
    }
  else{
    j0=j;
    j1=j+1;
    xi0=xi;
    xi1=grid.reduced_xi[j1];
    status=ind1D(&grid.x[xi1],grid.reduced_nx[j1],lon,&i1);
    }
  
  int i00,i01,i10,i11;
  
  if(lon<grid.x[xi0+i0]){
    i00=i0-1;
    i01=i0;
    }
  else{
    i00=i0;
    i01=i0+1;
    }
  
  if(lon<grid.x[xi1+i1]){
    i10=i1-1;
    i11=i1;
    }
  else{
    i10=i1;
    i11=i1+1;
    }
  
  /* wrap around */
  const int
    nx0=grid.reduced_nx[j0],
    nx1=grid.reduced_nx[j1];
  
  if(i00<0)i00=nx0-1;
  if(i01>=nx0)i01=0;
  if(i10<0)i10=nx1-1;
  if(i11>=nx1)i11=0;
  
  const int
    k00=xi0+i00,
    k01=xi0+i01,
    k10=xi1+i10,
    k11=xi1+i11;
  const double
    x00=degree_recale(grid.x[k00],lon),
    x01=degree_recale(grid.x[k01],lon),
    x10=degree_recale(grid.x[k10],lon),
    x11=degree_recale(grid.x[k11],lon),
    y0=grid.y[j0],
    y1=grid.y[j1];
  double wt=0.,w=0.,wm,wy;
  T z=0.,zm;
  
  wy=fabs(lat-y1);
  
  zm=buf[k00];
  wm=wy*fabs(lon-x01);
  incwz(zm,mask,wm,&wt,&w,&z);
  
  zm=buf[k01];
  wm=wy*fabs(lon-x00);
  incwz(zm,mask,wm,&wt,&w,&z);
  
  wy=fabs(lat-y0);
  
  zm=buf[k10];
  wm=wy*fabs(lon-x11);
  incwz(zm,mask,wm,&wt,&w,&z);
  
  zm=buf[k11];
  wm=wy*fabs(lon-x10);
  incwz(zm,mask,wm,&wt,&w,&z);
  
  if(w<wt*.5)
    TRAP_ERR_RETURN(mask,verbose-1,"close to masked %g<.5*%g\n",w,wt);
  
  z/=w;
  if(verbose>0) {
    STDERR_BASE_LINE_FUNC("");
    cerr << z << "(" << w << ")\n";
    }
  
  return z;
}


bool doBilinear=true;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void index_interpolation01(const grid_t &grid,double lon,double lat,int64_t *m,const T *buf,const T & mask,T *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///
/**
\param lon longitude in degrees
\param lat latitude in degrees
\param[in,out] *m index of the closest node. Used as prior solution.
*/
/*----------------------------------------------------------------------------*/
{
  int i,j,status;
  size_t k,xi;
  
  if(grid.modeH!=1 and grid.modeH!=MODEH_REDUCED_GG)
    TRAP_ERR_EXIT(ENOEXEC,"%s should not have been called with grid.modeH=%d\n",__func__,grid.modeH);
  
/*---------------------------------------------------------------------*//**<h1>
  dichotomy search of closest node </h1>*/
  
  if(*m!=-1 && *m<grid.Hsize()){
    grid.ij(*m,&i,&j);
    }
  else{
    i=-1;
    j=-1;
    }
  
  status=ind1D(grid.y,grid.ny,lat,&j);
  
  if(status!=0){
    *m=-1;
    if(z!=0)
      *z=mask;
    return;
    }
  
  if(grid.modeH==1){
    lon=map_recale(grid,lon);
    status=ind1D(grid.x,grid.nx,lon,&i);
    }
  else if(grid.modeH==MODEH_REDUCED_GG){
    const int nx=grid.reduced_nx[j];
    xi=grid.reduced_xi[j];
    degree_recale(&lon,(grid.x[xi]+grid.x[xi+nx-1])*.5);
    status=ind1D(&grid.x[xi],nx,lon,&i);
    status=0;
    }
  
  if(status!=0){
    *m=-1;
    if(z!=0)
      *z=mask;
    return;
    }
  
  if(grid.modeH==1){
    k=(size_t) j*(size_t) grid.nx+i;
    }
  else if(grid.modeH==MODEH_REDUCED_GG){
    k=xi+i;
    if(z!=0){
      if(doBilinear)
        *z=bilinear_reduced_gg_interpolation(grid,lon,lat,xi,i,j,buf,mask,verbose);
      else
        *z=buf[k];
      }
    }
  
  *m=k;
  if(z!=0 and grid.modeH!=MODEH_REDUCED_GG) {
    *z=buf[m4buf(grid,k)]; }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> int index_interpolation00(const grid_t &grid,double lon,double lat,const T *buf,const T & mask,T *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///
/*----------------------------------------------------------------------------*/
{
  int i,j;
  size_t k;
  
  if(grid.modeH!=0)TRAP_ERR_EXIT(ENOEXEC,"%s should not have been called with grid.modeH=%d\n",__func__,grid.modeH);
  
  i=round((lon-grid.xmin)/grid.dx);
  if(i<0 || i>=grid.nx) goto error;
  
  j=round((lat-grid.ymin)/grid.dy);
  if(j<0 || j>=grid.ny) goto error;
  
  k=j*grid.nx+i;
  if(z!=0)
    *z=buf[m4buf(grid,k)];
  return k;
  
error:
  if(z!=0)
    *z=mask;
  return -1;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int64_t get_neighbour(const grid_t &grid,int i,int j,int64_t m,int n)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///neighbour index
/**
\param grid
\param n 0 to 3 : right, top, left, bottom
\return index of neighbour or -1 if out of \c grid range
*/
/*----------------------------------------------------------------------------*/
{
  int64_t mn;
  
  while(n>=4) n-=4;
  while(n<0)  n+=4;
  
  switch(n){
  case 0: /* right */
    if(i>=grid.nx-1) return -1;
    mn=m+1;
    break;
  case 1: /* top */
    if(j==0) return -1;
    mn=m-grid.nx;
    break;
  case 2: /* left */
    if(i==0) return -1;
    mn=m-1;
    break;
  case 3: /* bottom */
    if(j>=grid.ny-1) return -1;
    mn=m+grid.nx;
    break;
  default:
    TRAP_ERR_EXIT(ENOEXEC,"coding error\n");
    }
  
  return mn;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> void incwz(T zm,T mask,double wm,double *wt,double *w,T *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  (*wt)+=wm;
  if(zm!=mask){
    (*w)+=wm;
    (*z)+=(T)wm*zm; /* THERE ARE NO operator*(float,complex<double>) */
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> T bilinear_interpolation(const grid_t &grid,double lon,double lat, int64_t m,const T *buf,const T & mask,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///bilinear part of interpolation
/**
\param lon longitude in degrees
\param lat latitude in degrees
\param m index of the closest node
*/
/*----------------------------------------------------------------------------*/
{
  if(m<0) TRAP_ERR_RETURN(mask,verbose-1,"wrong index\n");
  
  if(grid.modeH==MODEH_REDUCED_GG)
    TRAP_ERR_EXIT(ENOEXEC,"not coded for MODEH_REDUCED_GG\n");
  
  int i,j;
  vector2D_t v,v0(lon,lat),vm;
  
  i=m%grid.nx;
  j=m/grid.nx;
  grid.xy(m,vm.x,vm.y);
  
  v0.x=degree_recale(v0.x,vm.x);
  v=v0-vm;
  
  struct{
    vector2D_t v;
    int64_t m;
    int n;
    } v1,v2;
  double d1,d2,d12;
  bool p2;
  
  /* find quadrangle */
  for(v2.n=0;v2.n<=4;v2.n++){
    
    if(v2.n>0){
      v1=v2;
      d1=-d2;
      }
    
    /* unitary vector */
    v2.m=get_neighbour(grid,i,j,m,v2.n);
    if(v2.m<0) continue;
    grid.xy(v2.m,v2.v.x,v2.v.y);
    v2.v.x=degree_recale(v2.v.x,vm.x);
    v2.v-=vm;
    v2.v*=1./hypot(v2.v);
    d2=v*v2.v;
    if(v2.n<=0 || v1.m<0) continue;
    
    /* check whether between consecutive neighbours of closest node */
    /* i.e. d1, d2 and d12 have the same sign or 0 */
    if(!d2) break;
    p2=d2>0.;
    if(p2==(d1<0.)) continue;
    d12=v1.v*v2.v;
    if(p2==(d12<0.)) continue;
    
    break;
    }
  
  if(v2.n>4) TRAP_ERR_RETURN(mask,verbose-1,"outside: %d>4\n",v2.n);
  
  /* get node opposite to closest node in quadrangle */
  int64_t m_;
  vector2D_t v_,v3,v4;
  double d3,d4;
  
  i=v1.m%grid.nx;
  j=v1.m/grid.nx;
  m_=get_neighbour(grid,i,j,v1.m,v2.n);
  grid.xy(m_,v_.x,v_.y);
  v_.x=degree_recale(v_.x,vm.x);
  v=v0-v_;
  grid.xy(v2.m,v3.x,v3.y);
  v3.x=degree_recale(v3.x,vm.x);
  v3-=v_;
  v3*=1./hypot(v3);
  d3=v3*v;
  if(d2 && p2==(d3<0.)) TRAP_ERR_RETURN(mask,verbose-1,"outside: %g<>%g\n",d2,d3);
  grid.xy(v1.m,v4.x,v4.y);
  v4.x=degree_recale(v4.x,vm.x);
  v4-=v_;
  v4*=1./hypot(v4);
  d4=v*v4;
  if(d2 && p2==(d4<0.)) TRAP_ERR_RETURN(mask,verbose-1,"outside: %g<>%g\n",d2,d4);
  
  T z=0.,zm;
  double w=0.,wt=0.,wm;
  
  zm=buf[m4buf(grid,m)];
  if(verbose>0) {
    STDERR_BASE_LINE_FUNC("");
    cerr << zm << "\n";
    }
  wm=d3*d4;
  incwz(zm,mask,wm,&wt,&w,&z);
  zm=buf[m4buf(grid,v1.m)];
  wm=d3*d2;
  incwz(zm,mask,wm,&wt,&w,&z);
  zm=buf[m4buf(grid,v2.m)];
  wm=d1*d4;
  incwz(zm,mask,wm,&wt,&w,&z);
  zm=buf[m4buf(grid,m_)];
  wm=d1*d2;
  incwz(zm,mask,wm,&wt,&w,&z);
  
  if(w<wt*.5) TRAP_ERR_RETURN(mask,verbose-1,"close to masked %g<.5*%g\n",w,wt);
  
  z/=w;
  if(verbose>0) {
    STDERR_BASE_LINE_FUNC("");
    cerr << z << "(" << w << ")\n";
    }
  
  return z;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> void index_interpolation_template(const grid_t &grid,double lon,double lat, int64_t *m,const T *buf,const T & mask,T *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///wraper for index_interpolation00(), index_interpolation01() and index_interpolation02()
/**
\param lon longitude in degrees
\param lat latitude in degrees
\param[in,out] *m
*/
/*----------------------------------------------------------------------------*/
{
  switch(grid.modeH){
    case 0:
      *m=index_interpolation00(grid,lon,lat,buf,mask,z,verbose);
      break;
    case 1:
    case MODEH_REDUCED_GG:
      index_interpolation01(grid,lon,lat,m,buf,mask,z,verbose);
      break;
    case 2:
    case -2:
      if(grid.nx<=1 or grid.ny<=1){
        *z=mask;
        *m=-1;
        return;
        }
      if(grid.list==0) TRAP_ERR_EXIT(ENOEXEC,"%s should not have been called with grid_t::modeH==%d and grid_t::list=%p, i.e. without calling set_grid_list(grid_t*) first.\n",__func__,grid.modeH,grid.list);
      *m=index_interpolation02(grid,lon,lat,*grid.list,buf,mask,z,verbose);
      /* bilinear_interpolation() already called */
      return;
    default:
      TRAP_ERR_EXIT(ENOEXEC,"%s should not have been called with grid.modeH=%d\n",__func__,grid.modeH);
    }
  
  if(doBilinear and grid.nx>1 and grid.ny>1 and z!=0 and grid.modeH!=MODEH_REDUCED_GG)
    *z=bilinear_interpolation(grid,lon,lat,*m,buf,mask,verbose);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void index_interpolation(const grid_t & grid,double lon,double lat,int64_t *m,const float *buf,const float & mask,float *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  index_interpolation_template(grid,lon,lat,m,buf,mask,z,verbose);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void index_interpolation(const grid_t &grid,double lon,double lat,int64_t *m,const double *buf,const double & mask,double *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  index_interpolation_template(grid,lon,lat,m,buf,mask,z,verbose);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void index_interpolation(const grid_t &grid,double lon,double lat,int64_t *m,const double *buf,const double & mask,complex<double> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double d;
  
  index_interpolation_template(grid,lon,lat,m,buf,mask,&d,verbose);
  
  *z=d;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void index_interpolation(const grid_t &grid,double lon,double lat,int64_t *m,const complex<double> *buf,const complex<double> & mask,complex<double> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  index_interpolation_template(grid,lon,lat,m,buf,mask,z,verbose);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

void index_interpolation(const grid_t &grid,double lon,double lat,int64_t *m,const complex<float> *buf,const complex<float> & mask,complex<float> *z,int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  index_interpolation_template(grid,lon,lat,m,buf,mask,z,verbose);
}
