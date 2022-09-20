
/*******************************************************************************

  T-UGO tools, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France (PhD)
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

E-mail: florent.lyard@legos.obs-mip.fr

*******************************************************************************/

#include <config.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "tides.h"

#include "maths.h"
#include "rutin.h"    /*  rutin.h contains common utility routines  */
#include "geo.h"
#include "map.h"
#include "netcdf-proto.h"
#include "bmg.h"
#include "topo.h"
#include "statistic.h"


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  hconstant_t * load_atlas(const string & atlas_template,const string & v1,const string & v2, const spectrum_t & WaveList, int *pn_, poc_var_t *bv)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  hconstant_t *constants;
  
  poc_data_t<double> a,g;
  char *aP,*aN,*gN;/* decoded names of atlas,amplitudes and phases */
  
  int wI,status,pI,pn;
  
  for(wI=0;wI<WaveList.n;wI++){
    const tidal_wave *wave=&WaveList.waves[wI];
    
    tide_decode_atlasname(NULL,atlas_template.c_str(),wave->name, 0, &aP,3);
    tide_decode_atlasname(NULL,v1.c_str(),wave->name, 0, &aN,2);
    tide_decode_atlasname(NULL,v2.c_str(),wave->name, 0, &gN,2);
    /*GA=greenwhich_argument(astro_angles,WaveList.waves[i]);
    if(abs(GA<1e-99))GA=0;*/
    printf("Reading constants for %10s (%8.5g deg/h) from %s...",wave->name,wave->omega,aP);fflush(stdout);
    
    if(wI==0){
      status=a.init(aP,aN);
      status=g.init(aP,gN);
      if(bv!=0)
        bv->dimensions=a.info.dimensions;
      pn=a.length;
      if(pn_!=0)
        *pn_=pn;
      constants=new hconstant_t[pn];
      for(pI=0;pI<pn;pI++){
        constants[pI].init_polar(WaveList.n);
        }
      }
    else{
      a.info.name=aN;
      g.info.name=gN;
      }
    status=a.read_data(aP,0,0,1);
    status=g.read_data(aP,0,0,1);
    
    for(pI=0;pI<pn;pI++){
      double *aI=&a.data[pI];
      if(*aI==a.mask){
        *aI=NAN;
        }
      constants[pI].a[wI]=*aI;
      constants[pI].G[wI]=g.data[pI];
      }
    
    delete[]aP;
    delete[]aN;
    delete[]gN;
    
    printf("%d done\n",pn);
    }
  
  a.destroy_data();
  g.destroy_data();
  
  return constants;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int load_atlas(const char *filename,const char *v1,const char *v2, grid_t *grid, fcomplex **tide, fcomplex *cmask, int mode, int verbose)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/*----------------------------------------------------------------------------*/
{
  int i, j, n, status;
  float *buf[2],spec[2];
  fcomplex zz;
  cdfgbl_t global;
  int id;
  cdfvar_t info;

  status= cdf_globalinfo(filename,&global,0);
  if(status!=NC_NOERR){
    nc_check_error(status,__LINE__,__FILE__,"cdf_globalinfo error with %s",filename);
    return status;
    }
  id=cdf_identify(global,v1);
  if(id<0){
    status=NC_ENOTVAR;/* Variable not found */
    nc_check_error(status,__LINE__,__FILE__,"cdf_identify error with %s in %s",v1,filename);
    return status;
    }

  status= cdf_varinfo(filename,v1,&info,0);
  status=poc_getgrid2d  (filename, global, info, grid);
  if(status!=NC_NOERR){
    nc_check_error(status,__LINE__,__FILE__,"cdf_loadvargrid_2d error with %s in %s",v1,filename);
    return status;
    }

  exitIfNull(buf[0]  =new float   [grid->nx*grid->ny]);
  exitIfNull(buf[1]  =new float   [grid->nx*grid->ny]);
  exitIfNull(*tide   =new fcomplex[grid->nx*grid->ny]);

/*------------------------------------------------------------------------------
  load netcdf variable */
  id=cdf_identify(global,v1);
  status= cdf_varinfo(filename,v1,&info,0);
  status= poc_getvar2d (filename, id, 0, buf[0], &spec[0], info);
  if(status!=NC_NOERR){
    nc_check_error(status,__LINE__,__FILE__,"poc_getvar2d error with %s in %s",v1,filename);
    return status;
    }

  id=cdf_identify(global,v2);
  if(id<0){
    status=NC_ENOTVAR;/* Variable not found */
    nc_check_error(status,__LINE__,__FILE__,"error with %s in %s",v2,filename);
    return status;
    }
  status= cdf_varinfo(filename,v2,&info,0);
  status= poc_getvar2d (filename, id, 0, buf[1], &spec[1], info);
  if(status!=NC_NOERR){
    nc_check_error(status,__LINE__,__FILE__,"poc_getvar2d error with %s in %s",v1,filename);
    return status;
    }
  global.destroy();

  switch(mode) {
    case POLAR:
      *cmask=complex<float> (spec[0],spec[0]);
      if(isnan(spec[0])!=0) *cmask=complex<float>(1.e+10,1e+10);
      for (j=0;j<grid->ny;j++)
        for (i=0;i<grid->nx;i++) {
          n=i+grid->nx*j;
          if((buf[0][n]!=spec[0]) && (buf[1][n]!=spec[1])) {
            (*tide)[n]=polar<float>(buf[0][n],-buf[1][n]*d2r);
            }
         else {
           (*tide)[n]=*cmask;
           }
         }
       break;
    case CARTESIAN:
      *cmask=complex<float> (spec[0],spec[1]);
      if(isnan(spec[0])!=0) *cmask=complex<float>(1.e+10,1e+10);
      for (j=0;j<grid->ny;j++)
        for (i=0;i<grid->nx;i++) {
          n=i+grid->nx*j;
          if((buf[0][n]!=spec[0]) && (buf[1][n]!=spec[1]) && (isnan(buf[0][n])==0)) {
            (*tide)[n]=complex<float>(buf[0][n],-buf[1][n]);
            }
         else {
           (*tide)[n]=*cmask;
           }
         }
       break;
    }
  delete[]buf[0];
  delete[]buf[1];
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int swap_XY2YX_template(grid_t grid, T *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,m,n;
  T *tmp;

  tmp=new T[grid.ny*grid.nx];
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      n=i*grid.ny+j;
      tmp[m]=buffer[n];
      }
    }
  for(j=0;j<grid.ny;j++) {
    for(i=0;i<grid.nx;i++) {
      m=j*grid.nx+i;
      buffer[m]=tmp[m];
      }
    }
  delete[] tmp;
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int swap_XY2YX(grid_t grid, fcomplex *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=swap_XY2YX_template(grid, buffer);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int swap_XY2YX(grid_t grid, float *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=swap_XY2YX_template(grid, buffer);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int swap_XY2YX(grid_t grid, double *buffer)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;

  status=swap_XY2YX_template(grid, buffer);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int eot_loadc1(const char* filename,grid_t *grid, fcomplex **tide, fcomplex *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status;
  const char *v1="re",*v2="im";
  status=load_atlas(filename, v1, v2, grid, tide, mask,CARTESIAN,0);
  
  }
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int sirocco_loadc1(const char* filename, grid_t *grid, fcomplex **tide, fcomplex *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int status;
  const char *v1="sossheig_a",*v2="sossheig_G";
  status=load_atlas(filename, v1, v2, grid, tide, mask,POLAR,0);
  return(status);
  }
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int boy_loadc1(const char* filename, grid_t *grid, fcomplex **tide, fcomplex *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,j,m,nitems,status;
  double x,y;
  float a,b;
  FILE *in;
  
  *mask=complex<float>(9999.,9999.);
  
  
  grid->xmin=  0.;
  grid->xmax=360.;
  
  grid->ymin= -90+1./16.;
  grid->ymax= +90.;
  
  grid->dx=1./16.;
  grid->dy=1./16.;
  
  map_set2Dgrid(grid, 0., -90+1./16., 360., 90.0, 1./16., 1./16.);
  
  status=map_completegridaxis(grid,1);
  
  *tide=new complex<float> [grid->Hsize()];
  
  for(j=0; j<grid->ny; j++) {
    for(i=0; i< grid->nx; i++) {
      m=grid->nx*j+i;
      (*tide)[m]=*mask;
      }
    }
  
  in=fopen(filename,"r");
  
  for(j=grid->ny-1; j>0 ; j--) {
    for(i=0; i< grid->nx; i++) {
      nitems=fscanf(in,"%lf %lf %f %f",&x,&y,&a,&b);
      if(nitems!=4) {
        fclose(in);
        return(-1);
        }
      m=grid->nx*j+i;
      (*tide)[m]=polar<float>(a,-b*M_PI/180.0);
      }
    }
  
  fclose(in);
  
  return(status);
  }
  
  
//*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int osu_loadc1(const char* filename, grid_t *grid, fcomplex **tide, fcomplex *mask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,j,m,nitems,status;
  double x,y;
  float a,b;
  FILE *in;
  
  *mask=complex<float>(999.,-999.);
  
  
  grid->xmin=  0.+1./8.;
  grid->xmax=360.+1./8.;
  
  grid->ymin= -90+1./8.;
  grid->ymax= +90.-1./8.;
  
  grid->dx=1./4.;
  grid->dy=1./4.;
  
  map_set2Dgrid(grid, 0.+1./8., -90+1./8., 360.+1./8., 90.0-1./8., 1./4., 1./4.);
  
  status=map_completegridaxis(grid,1);
  
  *tide=new complex<float> [grid->Hsize()];
  
  for(j=0; j<grid->ny; j++) {
    for(i=0; i< grid->nx; i++) {
      m=grid->nx*j+i;
      (*tide)[m]=*mask;
      }
    }
  
  in=fopen(filename,"r");
  
  for(j=0; j<grid->ny; j++) {
    for(i=0; i< grid->nx-1; i++) {
      nitems=fscanf(in,"%lf %lf %f %f",&x,&y,&a,&b);
      if(nitems!=4) {
        fclose(in);
        return(-1);
        }
      m=grid->nx*j+i;
//      (*tide)[m]=polar<float>(a,-b*M_PI/180.0);
      (*tide)[m]=complex<float>(a,-b);
      }
    }
  
  for(j=0; j<grid->ny; j++) {
    m=grid->nx*j;
    (*tide)[m+grid->nx-1]=(*tide)[m];
    }
  
  fclose(in);
  
  return(status);
  }
  

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tpxo_loadc1(const char* filename, grid_t *grid,fcomplex ***tide, fcomplex *mask, int *nbuffers, spectrum_t *s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
  {
  int i,k,m,status;
  int id,ncid,dimid;
  int nwaves,cdim;
  float *realp, *imagp;
  cdfgbl_t global;
  cdfvar_t info;
  decoded_t decoded;
  range_t<double> range;
  char *names;
  
  status=cdf_globalinfo(filename,&global,0);
  
  dimid=cdf_identify_dimension(global,"nx");
  grid->nx=global.dimension[dimid].length;
  dimid=cdf_identify_dimension(global,"ny");
  grid->ny=global.dimension[dimid].length;
  
  grid->x=new double[grid->nx*grid->ny];
  grid->y=new double[grid->nx*grid->ny];
  
  status=nc_open(filename,0,&ncid);
  nc_check_error(status,__LINE__,__FILE__);
 
  status=cdf_varinfo(global,"lon_z",&info);

  switch(info.ndim) {
    case 1:
      grid->x=new double[grid->nx];
      status=nc_get_var_double(ncid,info.id,grid->x);
      nc_check_error(status,__LINE__,__FILE__);
      break;
    case 2:
      grid->x=new double[grid->nx*grid->ny];
      status=nc_get_var_double(ncid,info.id,grid->x);
      nc_check_error(status,__LINE__,__FILE__);
      status=swap_XY2YX(*grid, grid->x);
      break;
    }
//   id=cdf_identify(global,"lon_z");
//   status=nc_get_var_double(ncid,id,grid->x);
//   nc_check_error(status,__LINE__,__FILE__);
  
//   status=swap_XY2YX(*grid, grid->x);
  
  status=cdf_varinfo(global,"lat_z",&info);

  switch(info.ndim) {
    case 1:
      grid->y=new double[grid->ny];
      status=nc_get_var_double(ncid,info.id,grid->y);
      nc_check_error(status,__LINE__,__FILE__);
      grid->modeH=1;
      break;
    case 2:
      grid->y=new double[grid->nx*grid->ny];
      status=nc_get_var_double(ncid,info.id,grid->y);
      nc_check_error(status,__LINE__,__FILE__);
      status=swap_XY2YX(*grid, grid->y);
      grid->modeH=2;
      break;
    }
  
//   id=cdf_identify(global,"lat_z");
//   status=nc_get_var_double(ncid,id,grid->y);
//   nc_check_error(status,__LINE__,__FILE__);
  
//   status=swap_XY2YX(*grid, grid->y);
  
  map_minmax(grid);
  
  grid->dx=(grid->xmax-grid->xmin)/((double) grid->nx-1.);
  grid->dy=(grid->ymax-grid->ymin)/((double) grid->ny-1.);

//   grid->modeH=2;
  
  
  dimid=cdf_identify_dimension(global,"nc");
  
  if(dimid==-1) {
    nwaves=1;
    }
  else {
    nwaves=global.dimension[dimid].length;
    }
  dimid=cdf_identify_dimension(global,"nct");
  cdim=global.dimension[dimid].length;
  
  status=cdf_varinfo(global,"con",&info);
  
  id=cdf_identify(global,"con");
  
  names=new char[nwaves*cdim];
  status=nc_get_var_text(ncid,id,names);
  nc_check_error(status,__LINE__,__FILE__);
  info.destroy();
  
  
  realp=new float[grid->nx*grid->ny*nwaves];
  imagp=new float[grid->nx*grid->ny*nwaves];
   
  status=cdf_varinfo(global,"hRe",&info);
  
  status=nc_get_var_float(ncid,info.id,realp);
  nc_check_error(status,__LINE__,__FILE__);
  info.destroy();
   
  status=cdf_varinfo(global,"hIm",&info);
  
  status=nc_get_var_float(ncid,info.id,imagp);
  nc_check_error(status,__LINE__,__FILE__);
  info.destroy();
  
  *tide=new fcomplex *[nwaves];
  int size=grid->nx*grid->ny;
  for(k=0;k<nwaves;k++) {
    (*tide)[k]=new fcomplex[grid->nx*grid->ny];
    for(m=0;m<grid->nx*grid->ny;m++) {
      (*tide)[k][m]=complex<float>(realp[k*size+m],imagp[k*size+m]);
      }
    status=swap_XY2YX(*grid, (*tide)[k]);
    }
  status = nc_close(ncid);
  nc_check_error(status,__LINE__,__FILE__);

  *mask=complex<float>(0.,0.);
  
  *nbuffers=nwaves;
  
  s->n=nwaves;
  s->waves=new tidal_wave[s->n];
  for(k=0;k<nwaves;k++) {
    char tmp[4];
    strncpy(tmp,(const char*) &(names[k*cdim]),4);
    for(i=0;i<cdim;i++) if(tmp[i]==' ') tmp[i]=0;
    strcpy(s->waves[k].name,tmp);
    }
  return(0);
  }



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ellipse_Madmp(const complex<double> u,const complex<double> v,double *M,double *a,double *d,double *m,double *p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///convert UV currents to ellipses
/**
\param u
\param v
\param *M maximum current
\param *a phase of maximum current in degrees
\param *d direction of maximum current in degrees
\param *m minimum current
\param *p polarisation: set to 1. if clockwise
See also ellipse_M()
*/
/*----------------------------------------------------------------------------*/
{
  const complex<double>
    u2=u*u,
    v2=v*v,
    u2pv2=u2+v2;
  const double
    au2pav2=abs(u2)+abs(v2),
    au2pv2=abs(u2pv2);
  
  if(M!=0)
    *M=sqrt((au2pav2+au2pv2)*.5);
  
  if(a!=0 || d!=0){
    const double ar=arg(u2pv2)*.5;
    if(d!=0){
      *d=atan2(
          abs(u)*cos(arg(u)-ar),
          abs(v)*cos(arg(v)-ar)
        )*r2d;
      if(*d<0.)*d+=360.;
      }
    if(a!=0)
      *a=-ar*r2d;
    }
  
  if(m!=0)
    *m=sqrt((au2pav2-au2pv2)*.5);
  
  if(p!=0)
    /* vectorial product */
    *p=sign(imag(u)*real(v)-real(u)*imag(v));
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ellipse_Madmp3D(const complex<double> u,const complex<double> v,const complex<double> w,double *M,double *a)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// convert UVW currents to ellipses
/**
\param u
\param v
\param w
\param *M maximum current
See also ellipse_M()
*/
/*----------------------------------------------------------------------------*/
{
  const complex<double>
    u2=u*u,
    v2=v*v,
    w2=w*w,
    sumS=u2+v2+w2;
  const double
    sumA=abs(u2)+abs(v2)+abs(w2),
    aSumS=abs(sumS);
  
  if(M!=0)
    *M=sqrt((sumA+aSumS)*.5);
  
  if(a!=0){
    const double ar=arg(sumS)*.5;
    if(a!=0)
      *a=-ar*r2d;
    }
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  template<typename T> ellipse_t ellipse_parameter_template(const complex<T> vx,const complex<T> vy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  complex<double> z;
  double x,y;
  double p1,p2,q1,q2,l,n,lpn,lmn,a,m,w,pol;
  ellipse_t e;
  
/*------------------------------------------------------------------------------
     Rmin: minor axis of the vector ellipse
     Rmax: major axis of the vector ellipse
     Pol : polarisation of the ellipse = rotation of the vector
           1 if direct, -1 otherwise
     Dir : inclination of the major axis on the WE axis
     Time: phase lag between phase origin and phase of maximum amplitude

     p=wt+V0

     a=Ax cos(-Gx), b=Ax sin(-Gx)

     X=Ax cos(p-Gx)= a cos(p) - b sin(p)
     Y=Ay cos(p-Gy)= c cos(p) - d sin(p)
     
     L=a²+c²
     N=b²+d²
     M=ab+cd

           Ro²=(X²+Y²) = (a²+c²) cs² -2*(ab+cd) cs*sn +(b²+d²2) sn²
                       = L*cs²-2*M*cs*sn+N*sn²
                       = (L-2*Mt+Nt²)/(1+t²)
                       = 0.5*[(L+N)-2*Msin2p+(L-N)cos2p]

     Ro extrema

     dRo/dp=0 -> 2*(N-L)*cs*sn-2*M*(cs²-sn²)=0

                2*cs*sn    =sin(2p)
                cs²-sn²    =cos(2p)

     dRo/dp=0 -> (N-L)*sin(2p)-2*M*cos(2p)=0
     
                 M*t²+(N-L)*t-M=0

                 t= -2K + +/-sqrt(K²+4)   K=(N-L)/M

                 t1={-(N-L)-sqrt[(N-L)²+4*M²]}/2M
                 t2={-(N-L)+sqrt[(N-L)²+4*M²]}/2M

     if N-L#0 -> tg(2p)=2M/(N-L)
                 cs=1-t²/1+t²

     if N-L=0 : if M#0 -> p=+/- 90.0
              : if M=0 -> p=undefined

-------------------------------------------------------------------------------

           Ro²=(X²+Y²)= 0.5*[(L+N)-2Msin2p+(L-N)cos2p]
                      = 0.5*[(L+N)+C cos(2p-W)]

     if L-N#0:

                     C=sqrt[4*M²+(L-N)²]

                Rmax² = 0.5*[(L+N)+C]    p=0.5*W     +/-pi
                Rmin² = 0.5*[(L+N)-C]    p=0.5*[W-pi]+/-pi
               cos(W) = (L-N)/C
               sin(W) = -2M/C
               tan(W) = -2M/(L-N)

     if L-N>0: -pi/2 < W < +pi/2
     if L-N<0: +pi/2 < W < +3*pi/2

     
     if L-N=0:


                   Ro²= 0.5*[(L+N)-2Msin2p]
                      = 0.5*[(L+N)-2Mcos(pi/2-2p)]

                     C=2*abs(M)

                Rmax² = 0.5*[(L+N)+C]
                Rmin² = 0.5*[(L+N)-C]
     if M>0:
                    W = -pi/2

     if M=0:
                    W undefined

     if M<0:
                    W = +pi/2

     Inclination:
                    Vmax=[acos(0.5*[W-pi])-bsin(0.5*[W-pi]),
                          ccos(0.5*[W-pi])-dsin(0.5*[W-pi])]

     Polarisation:


!-----------------------------------------------------------------------------*/

  p1=real(vx);
  p2=imag(vx);

  q1=real(vy);
  q2=imag(vy);

  l=p1*p1+q1*q1;
  n=p2*p2+q2*q2;
  m=p1*p2+q1*q2;

  lmn=l-n;
  lpn=l+n;

  a=sqrt(lmn*lmn+4.*m*m);

  e.a=sqrt(0.5*fabs(lpn+a));
  e.b=sqrt(0.5*fabs(lpn-a));

  pol=p2*q1-p1*q2;
  e.polarisation=sign(pol);

  if(a != 0.0) {
    w=atan2(-2*m,lmn);
    e.phase=0.5*w;
    }
  else {
    e.phase=0.0;
    }

  z=polar<double>(1,e.phase);

  x=real((complex<double>) vx*z);
  y=real((complex<double>) vy*z);

/*------------------------------------------------------------------------*//**
  L'inclinaison de l'ellipse est calculee par rapport a l'axe
  Nord/Sud, clockwise */
/*----------------------------------------------------------------------------*/
  e.inclination=M_PI_2-atan2(y,x);
  if(e.inclination<0.) e.inclination+=M_PIx2;

  return(e);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  ellipse_t ellipse_parameter(const complex<float> vx, const complex<float> vy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ellipse_t e;
  
  e=ellipse_parameter_template(vx, vy);
  
  return(e);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  ellipse_t ellipse_parameter(const complex<double> vx, const complex<double> vy)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ellipse_t e;
  
  e=ellipse_parameter_template(vx, vy);
  
  return(e);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  void ellipse_parameter(fcomplex vx,fcomplex vy,float *rmin,float *rmax,float *pol, float *dir, float *time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  ellipse_t e;
  
  *rmax=e.a;
  *rmin=e.b;
  
  *pol=e.polarisation;
  
  *time=e.phase;

/*------------------------------------------------------------------------*//**
  L'inclinaison de l'ellipse est calculee par rapport a l'axe
  Ouest->Est, anti-clockwise */
/*----------------------------------------------------------------------------*/
  *dir=-e.inclination-M_PI_2;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tides_savec1(const char *output, grid_t grid, fcomplex *tide, fcomplex cmask, float scale, string & Aname, string & Gname)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i,j,n,status;
  fcomplex z;
  cdfvar_t variable;
  pocgrd_t ncgrid;
  float *rbufx=NULL,*rbufy=NULL, rmask;

  printf("#################################################################\n");
  printf("convert complex to amplitude, phase lag, units scale=%f\n",scale);
  rmask=cmask.real();
  
  exitIfNull(
    rbufx=new float[grid.nx*grid.ny]
    );
  exitIfNull(
    rbufy=new float[grid.nx*grid.ny]
    );
    
  for (j=0;j<grid.ny;j++) {
    for (i=0;i<grid.nx;i++) {
      n=i+grid.nx*j;
      z=tide[n];
      if(z!=cmask){
        rbufx[n]=scale*abs(z);
        rbufy[n]=-arg(z)*r2d;
        }
      else {
        rbufx[n]=rmask;
        rbufy[n]=rmask;
        }
      }
    }
  
  printf("#################################################################\n");
  printf("write output file : %s\n",output);
  status= poc_createfile(output);
  grid.nz=1;
  grid.z=NULL;

  status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);

  poc_standardvariable_xy(&variable,Aname.c_str(),rmask,"m",1., 0.,"tidal_amplitude","tidal_amplitude","tidal_amplitude",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  grid, variable.id,rbufx);
  variable.destroy();

  poc_standardvariable_xy(&variable,Gname.c_str(),rmask,"degree",1., 0.,"tidal_phase_lag","tidal_phase_lag","tidal_phase_lag",ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output,  grid, variable.id,rbufy);
  variable.destroy();

//   if(topofile!=NULL) {
//     poc_standardvariable_xy(&variable, "bathymetry",1e+11,"m",1., 0.,"bathymetry","bathymetry","bathymetry",ncgrid);
//     status=create_ncvariable(output, &variable);
//     status=poc_write_xy(output,  topogrid, variable.id,topo);
//     variable.destroy();
//     }
  free(rbufx);
  free(rbufy);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tides_savec1(const char *output, grid_t grid, fcomplex *tide, fcomplex cmask, float scale)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  string Aname, Gname;
  
  Aname="Ha";
  Gname="Hg";
  
  status=tides_savec1(output, grid, tide, cmask, scale, Aname, Gname);
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tides_savec1(const char *output, const char *v1, const char *v2, const char *name, const char *units, grid_t zgrid, fcomplex *cbuf, fcomplex cmask, pocgrd_t ncgrid)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int i, status;
  float *rbufx, *rbufy;
  float rmask;
  fcomplex z;
  cdfvar_t variable;
  char standardname[1024], longname[1024];
  
  rmask=real(cmask);
  rbufx =new float[zgrid.nx*zgrid.ny];
  rbufy =new float[zgrid.nx*zgrid.ny];
  for (i=0;i<zgrid.nx*zgrid.ny;i++) {
    if(cbuf[i]!=cmask) {
      z=cbuf[i];
/*------------------------------------------------------------------------------
      convert fcomplex h in amplitude (meters) and phase lag (degrees)*/
      rbufx[i]=abs(z);
      rbufy[i]=-arg(z)*r2d;
      if(rbufy[i]<0.0) rbufy[i]+=360.0;
      }
    else {
      rbufx[i]=rmask;
      rbufy[i]=rmask;
      }
    }

  sprintf(standardname,"tidal_%s_amplitude",name);
  for(int k=0;k<strlen(standardname); k++) if(standardname[k]==' ') standardname[k]='_';
  sprintf(longname,"tidal %s amplitude",name);
  poc_standardvariable_xy(&variable,v1,rmask,units,1., 0., standardname,longname,v1,ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output, zgrid, variable.id,rbufx);
  variable.destroy();

  sprintf(standardname,"tidal_%s_phaselag",name);
  for(int k=0;k<strlen(standardname); k++) if(standardname[k]==' ') standardname[k]='_';
  sprintf(longname,"tidal %s phaselag",name);
  poc_standardvariable_xy(&variable,v2,rmask,"degree",1., 0.,standardname,longname,v2,ncgrid);
  status=create_ncvariable(output, &variable);
  status=poc_write_xy(output, zgrid, variable.id,rbufy);
  variable.destroy();

  delete[] rbufx;
  delete[] rbufy;
    
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tides_savec1(const char *output, const char *v1, const char *v2, const char *name, const char *units, grid_t grid, fcomplex *cbuf, fcomplex cmask)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  pocgrd_t ncgrid;
  
  status= poc_createfile(output);
  status=poc_sphericalgrid_xy(output,"",grid,&ncgrid);
  status=tides_savec1(output, v1, v2, name, units, grid, cbuf, cmask, ncgrid);
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int tides_load()
/* unused */

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  float  time;
  int nbuffers,status;
  char *input=NULL;
  char format[1024]="bmg";
  grid_t grid,atlas_grid, topogrid;
  fcomplex **tide=NULL,cmask,z,e;
  spectrum_t spectrum;
  bool extend;
  
  if(input == NULL) goto error;

  printf("#################################################################\n");
  printf("load input file : %s, format %s\n",input, format);
  if(strcmp(format,"bmg")==0) {
    tide=new fcomplex *[1];
    spectrum.n=1;
    status= bmg_loadc1(input,1,1,1,&atlas_grid,&(tide[0]),&cmask,&time);
    nbuffers=1;
    }
  else if(strcmp(format,"ascii")==0 or strcmp(format,"dtu")==0) {
    nbuffers=1;
    tide=new fcomplex *[nbuffers];
    spectrum.n=1;
    if(nbuffers==1)
      status= ascii_loadc2(input,&atlas_grid,&tide[0],0,&cmask);
    else
      status= ascii_loadc2(input,&atlas_grid,&tide[0],&tide[1],&cmask);
    }
  else if(strcmp(format,"netcdf")==0) {
    tide=new fcomplex *[1];
    spectrum.n=1;
    status= sirocco_loadc1(input,&atlas_grid,&(tide[0]),&cmask);
    nbuffers=1;
    }
  else if(strcmp(format,"got")==0) {
    tide=new fcomplex *[1];
    spectrum.n=1;
    status= ascii_loadc1_got(input,&atlas_grid,&(tide[0]),&cmask);
    nbuffers=1;
    }
  else if(strcmp(format,"tpxo")==0) {
    status=tpxo_loadc1(input,&atlas_grid,&tide,&cmask,&nbuffers,&spectrum);
//    spectrum.n=1;
    }
  else if(strcmp(format,"otis")==0) {
//    status= tpxo_loadc1(input,&atlas_grid,&tide,&cmask);
    }
  else if(strcmp(format,"eot")==0) {
    tide=new fcomplex *[1];
    spectrum.n=1;
    status= eot_loadc1(input,&atlas_grid,&(tide[0]),&cmask);
    nbuffers=1;
    }
  else {
    goto error;
    }

  if(status != 0) goto error;

  status=map_completegridaxis(&atlas_grid,2);

  if(extend) {
    int modified=0;
    if(nbuffers==1) {
      status=map_extendfield(&atlas_grid, &(tide[0]), &modified);
      }
    else {
      status=map_extendfield(&atlas_grid, tide, nbuffers, &modified);
      }
    }
  

  
  STDOUT_BASE_LINE("end of tides-converter... \n");
  exit(0);
error:
  STDOUT_BASE_LINE("error detected, quit ... \n");
  exit(-1);
}
