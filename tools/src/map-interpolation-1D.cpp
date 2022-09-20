
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief 1D interpolation definitions
*/
/*----------------------------------------------------------------------------*/

#include <complex>
#include <vector>
#include <algorithm> // for min and max

#include "fe.def" /* for LGP0 & co. */
#include "maths.h" //for vector3_t

using namespace std;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename A,typename T> int incrementi(const A & t,int n,int *i,int di,T time)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///step towards the right point
/**
\param di step size
\returns 0 when the right point is within the boundaries, whether a step has been taken or not, -1 or 1 when outside and closest to the lower or upper bounds
*/
/*----------------------------------------------------------------------------*/
{
  double l,lmin=+INFINITY,factor=1.;
  int ig,imin,il;
  
  if(di<0)di=-di;
  if(t[max(*i-di,0)]>t[min(*i+di,n-1)])factor=-1.;
  
  ig=max(*i-di,0);
  
  if(factor*time<factor*t[ig]){
    *i=ig;
    return -1;
    }
  
  ig=min(*i+di,n-1);
  
  if(factor*time>factor*t[ig]){
    *i=ig;
    return +1;
    }

  if(di==0)return 0;
  
  il=0;
  for(ig=*i-di;ig<n && il<3;ig+=di,il++){
    if(ig<0)continue;
    
    l=fabs(t[ig]-time);
    
    if(lmin>l){
      lmin=l;
      imin=ig;
      }
    
    }
  
  *i=imin;
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename A,typename T> int ind1D_template(const A & t,int n,T time,int *m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///find index of closest value
/**
\param t array[\c n] of sorted values
\param n
\param time
\param[in,out] *m index of the closest value. Used as prior solution.
\returns 0 when inside, -1 or 1 when outside
*/
/*----------------------------------------------------------------------------*/
{
  int i,di,status;
  const int n_2=n/2;
  
/*---------------------------------------------------------------------*//**<h1>
  dichotomy search of closest index </h1>*/
  
  if(*m>=0 && *m<n){
    di=1;
    i=*m;
    }
  else{
    redo:
    di=n_2;
    i=n_2;
    }
  
  //accelerate
  do{
    status=incrementi(t,n,&i,di,time);
    if( status==0 or di==n_2 or
        i<=0 or i>=n-1 )
      break;
    di*=2;
    if(di>n_2)
      goto redo;
    }while(1);
  
  if(di==n_2 && status!=0){
    //outside
    *m=i;
    return status;
    }
  
  //slow down
  while(di>1){
    di=(di+1)/2;
    status=incrementi(t,n,&i,di,time);
    }
  
  *m=i;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ind1D(const int *t,int n,int time,int *m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=ind1D_template(t,n,time,m);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ind1D(const float *t,int n,float time,int *m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=ind1D_template(t,n,time,m);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ind1D(const double *t,int n,double time,int *m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=ind1D_template(t,n,time,m);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int ind1D(const vector<double> & t,int n,double time,int *m)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int result;
  
  result=ind1D_template(t,n,time,m);
  
  return result;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T1, typename T2, typename V, typename M> int map_interpolate1D_template(const V & h,const T1 t, M mask, int n, T2 time, M *z, T2 tmax, u_int8_t extrapolate=0, int *k0=0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// 1D interpolation
/**
\param h array[\c n] of interpolated values
\param t array[\c n] of sorted values
\param mask
\param n
\param time at which to interpolate
\param *z result
\param tmax maximum interval over which to interpolate. If -1: disabled
\param extrapolate Ox10:nearest neighbour; 0x1:extrapolation; can also have discretisation (see code below)
\param[in,out] *kO index of the closest value. Used as prior solution for the search.
\returns almost always 0. -1 if a NAN value as been met in t
*/
/*----------------------------------------------------------------------------*/
{
  int k,k1,k2;
  T2 r;
  
  const u_int8_t mode=(extrapolate>>4)-1;
  extrapolate &= 0xf;

  if(time < t[0]) {
    if(extrapolate)
      *z=h[0];
    else
      *z=mask;
    return(0);
    }

  if(time > t[n-1]) {
    if(extrapolate)
      *z=h[n-1];
    else
      *z=mask;
    return(0);
    }

  if(k0!=0){
    ind1D(t,n,time,k0);
    k1=max(0,*k0-1);
    k2=min(*k0+1,n-1);
    }
  else{
    k1=0;
    k2=n-1;
    }

  while(k2-k1>1) {
    k=int( 0.5*(k2+k1) );
    if(time > t[k]) {
      k1=k;
      continue;
      }
    if(time < t[k]) {
      k2=k;
      continue;
      }
    if(time == t[k]) {
      *z=h[k];
      return(0);
      }
    *z=mask;
    return -1;
    }
  
  if(time == t[k2]) {
    *z=h[k2];
    return(0);
    }
  
  if(tmax !=-1) {
    if(t[k2]-t[k1]>tmax) {
      *z=mask;
      return(0);
      }
    }
  
  switch(mode){
  case 0:/* nearest neighbour or LGP0 */
    *z=h[k1];
    break;
  default:
    r=(t[k2]-time)/(t[k2]-t[k1]);
    
    if(mode==DGP1){
      k1*=2;
      k2=k1+1;
      }
    
    if(h[k1]==mask) {
      if(r>0.5)
        *z=mask;
      else
        *z=h[k2];
      return(0);
      }
    if(h[k2]==mask) {
      if(r<0.5)
        *z=mask;
      else
        *z=h[k1];
      return(0);
      }
    *z=r*h[k1]+(1.-r)*h[k2];
    }
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const float *h,const float *t,float mask,int n,float time,float *z,int extrapolate, int *k0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  float tmax=-1.0;
  
  status=map_interpolate1D_template(h, t, mask, n, time, z, tmax, extrapolate, k0);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const float *h,const double *t,float mask,int n,double time,float *z,int extrapolate, int *k0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  double tmax=-1.0;
  
  status=map_interpolate1D_template(h, t, mask, n, time, z, tmax, extrapolate, k0);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const vector<float> & h,const vector<double> & t, float mask, int n, double time, float *z, int extrapolate, int *k0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  double tmax=-1.0;
  
  status=map_interpolate1D_template(h, t, mask, n, time, z, tmax, extrapolate, k0);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const double *h,const double *t, double mask, int n, double time, double *z, int extrapolate, int *k0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  double tmax=-1.0;
  
  status=map_interpolate1D_template(h, t, mask, n, time, z, tmax, extrapolate, k0);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const double *h,const double *t,double mask,int n,double time,complex<double> *z,int extrapolate, int *k0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  double d,tmax=-1.0;
  
  status=map_interpolate1D_template(h, t, mask, n, time, &d, tmax, extrapolate, k0);
  
  *z=d;
  
  return status;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const complex<double> *h,const double *t,complex<double> mask,int n,double time,complex<double> *z,int extrapolate, int *k0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  double tmax=-1.0;
  
  status=map_interpolate1D_template(h, t, mask, n, time, z, tmax, extrapolate, k0);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const double *h,const double *t, double mask, int n, double time, double *z, double tmax, int extrapolate, int *k0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  
  status=map_interpolate1D_template(h, t, mask, n, time, z, tmax, extrapolate, k0);
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const double *h,const double *t, double mask, int n, double time, short unsigned int *z, double tmax, int extrapolate, int *k0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  double zz;
  status=map_interpolate1D_template(h, t, mask, n, time, &zz, tmax, extrapolate, k0);
  *z=zz;
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const vector3_t *h,const double *t, int n, double time, vector3_t *z, int *k0)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status=0;
  static const vector3_t mask(NAN,NAN,NAN);
  
  status=map_interpolate1D_template(h, t, mask, n, time, z, -1., 0, k0);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int map_interpolate1D(const complex<float> *h,const double *t, complex<float> mask, int n, double time, complex<float> *z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int k,k1,k2,status=0;
  float r=0;
  complex<float> c;

//   if(time < t[0]) {
//     *z=mask;
//     return(0);
//     }
//
//   if(time > t[n-1]) {
//     *z=mask;
//     return(0);
//     }

  if((time - t[0])*(time - t[n-1]) > 0.0) {
    *z=mask;
    return(0);
    }

  k1=0;
  k2=n-1;

  if(t[0]<t[n-1]) {
  while(k2-k1>1) {
    k=int( 0.5*(k2+k1) );
    if(time > t[k]) {
      k1=k;
      continue;
      }
    if(time < t[k]) {
      k2=k;
      continue;
      }
    if(time = t[k]) {
      *z=h[k];
      return(0);
      }
    }
    }
  else {
  while(k2-k1>1) {
    k=int( 0.5*(k2+k1) );
    if(time < t[k]) {
      k1=k;
      continue;
      }
    if(time > t[k]) {
      k2=k;
      continue;
      }
    if(time = t[k]) {
      *z=h[k];
      return(0);
      }
    }
    }

  if(h[k1]==mask) {
    *z=mask;
    return(0);
    }

  if(h[k2]==mask) {
    *z=mask;
    return(0);
    }

  r=(t[k2]-time)/(t[k2]-t[k1]);

  c=r*h[k1]+(1.f-r)*h[k2];

  *z=c;

  return(status);
}
