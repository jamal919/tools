
/*******************************************************************************

  T-UGO tools, 2006-2017

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\brief geometrical classes and function prototypes
*/
/*----------------------------------------------------------------------------*/

#ifndef MATHS_H
#define MATHS_H

#include <stdio.h>
#include <complex>
#include <vector>

#include "poc-assertions.h"
#include "errno.h"

using namespace std;


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> bool isfinite(const complex<T> & z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  const bool
    fr=isfinite(real(z)),
    fi=isfinite(imag(z)),
    f=fr && fi;
  
  return f;
}


template <typename T,int N>  class matrix_template_t {
private:
public:
  T c[N][N];
  
  matrix_template_t<T,N> () {
    int i,j;
    for(j=0;j<N;j++)
      for(i=0;i<N;i++)
        this->c[j][i]=0;
    }
  
  matrix_template_t<T,N> (T a, T b, T c, T d) {
    if(N!=2) TRAP_ERR_EXIT(ENOEXEC,"Coding error\n");
    /*  / a b \
        \ c d /  */
    this->c[0][0]=a;
    this->c[0][1]=b;
    this->c[1][0]=c;
    this->c[1][1]=d;
    }
  
  matrix_template_t<T,N> & operator*=(double k){
    int i,j;
    for(j=0;j<N;j++)
      for(i=0;i<N;i++)
        this->c[j][i]*=k;
    return *this;
    }
  
  matrix_template_t<T,N> operator * (const matrix_template_t<T,N> & __x) {
    matrix_template_t<T,N> __r;
    int i,j,k;
    
    for(j=0;j<N;j++)
      for(i=0;i<N;i++)
        for(k=0;k<N;k++)
          __r.c[j][i]+=c[j][k]*__x.c[k][i];
    
    return __r;
    }

  matrix_template_t<T,N> inverse() {
    matrix_template_t<T,N> __r;
    T det;
    if(N!=2) TRAP_ERR_EXIT(ENOEXEC,"Not coded yet for N!=2\n");
    det=c[0][0]*c[1][1]-c[0][1]*c[1][0];
    if(abs(det)==0.) {
      printf("Singluar matrix, SHOULD abort.\n");
      }
    __r.c[0][0]= c[1][1]/det;
    __r.c[0][1]=-c[0][1]/det;
    __r.c[1][0]=-c[1][0]/det;
    __r.c[1][1]= c[0][0]/det;
    return __r;
    }
};

typedef matrix_template_t<double,2> matrix2x2_t;
typedef matrix_template_t<complex<double>,2> zmatrix2x2_t;

typedef matrix_template_t<double,3> matrix3x3_t;
typedef matrix_template_t<complex<double>,3> zmatrix3x3_t;

typedef vector<int> cycle_t;

class point2D_t
{
private:
public:
  double  x,y;
  
  void init(double x0=0.,double y0=0.){
    x=x0;
    y=y0;
    }
  
  point2D_t(double x0=NAN,double y0=NAN){
    init(x0,y0);
    }
  
  inline point2D_t & operator+=(const point2D_t & a){
    x+=a.x;
    y+=a.y;
    return *this;
    }
  
  inline point2D_t & operator/=(const double a){
    x/=a;
    y/=a;
    return *this;
    }

  };

template<typename T0> class vector3_template_t {
public:
  T0 x,y,z;

  void init(){
    x=0.;
    y=0.;
    z=0.;
    }
  void init(T0 xx,T0 yy,T0 zz){
    x=xx;
    y=yy;
    z=zz;
    }
  
  template<typename T> inline vector3_template_t& operator *= (const T c) {
    x *=  c;
    y *=  c;
    z *=  c;
    return *this;
    }
  
  template<typename T> inline vector3_template_t& operator /= (const T c) {
    x /=  c;
    y /=  c;
    z /=  c;
    return *this;
    }
  
  inline vector3_template_t & operator-=(const vector3_template_t & a){
    x-=a.x;
    y-=a.y;
    z-=a.z;
    return *this;
    }
  
  inline vector3_template_t & operator+=(const vector3_template_t & a){
    x+=a.x;
    y+=a.y;
    z+=a.z;
    return *this;
    }
  
  };

class vector3_t : public vector3_template_t<double> {
public:

  vector3_t(){
    init();
    }
  vector3_t(double xx,double yy,double zz){
    init(xx,yy,zz);
    }
  
  };

  inline double square(const vector3_t & a){
    return a.x*a.x+a.y*a.y+a.z*a.z;
    }

  inline double manhattan(const vector3_t & a){
    return fabs(a.x)+fabs(a.y)+fabs(a.z);
    }

  inline double hypot(const vector3_t & a){
    return sqrt(square(a));
    }

  //* operator
  //scalar product
  inline double operator*(const vector3_t & a,const vector3_t & b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
    }

  inline vector3_t operator*(const vector3_t & a,const double k){
    return vector3_t(a.x*k,a.y*k,a.z*k);
    }

  inline vector3_t operator*(const double k,const vector3_t & a){
    return vector3_t(a.x*k,a.y*k,a.z*k);
    }

  //- operator
  inline vector3_t operator-(const vector3_t & a,const vector3_t & b){
    return vector3_t(a.x-b.x,a.y-b.y,a.z-b.z);
    }
  
  //+ operator
  inline vector3_t operator+(const vector3_t & a,const vector3_t & b){
    return vector3_t(a.x+b.x,a.y+b.y,a.z+b.z);
    }
  
  //& operator
  //vectorial product
  inline vector3_t operator&(const vector3_t & a,const vector3_t & b){
    return vector3_t(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
    }

  //^ operator
  //angle from 2 vectors
  inline double operator^(const vector3_t & a,const vector3_t & b){
    return atan2(hypot(a&b),a*b);
    }
  
  inline bool operator==(const vector3_t & a,const vector3_t & b){
    return a.x==b.x and a.y==b.y and a.z==b.z ;
    }

class cvector3_t : public vector3_template_t<complex<double> > {
public:

  cvector3_t(){
    init();
    }
  cvector3_t(complex<double> xx,complex<double> yy,complex<double> zz){
    init(xx,yy,zz);
    }
  
  inline cvector3_t & operator*=(double k){
    x*=k;
    y*=k;
    z*=k;
    return *this;
    }
  
  };
  
  //+ operator
  inline cvector3_t operator+(const cvector3_t & a,const cvector3_t & b){
    return cvector3_t(a.x+b.x,a.y+b.y,a.z+b.z);
    }

  //* operator
  //scalar product
  inline complex<double> operator*(const cvector3_t & a,const vector3_t & b){
    return a.x*b.x+a.y*b.y+a.z*b.z;
    }

  inline cvector3_t operator*(const complex<double> k,const vector3_t & a){
    return cvector3_t(a.x*k,a.y*k,a.z*k);
    }
  inline cvector3_t operator*(const vector3_t & a,const complex<double> k){
    return cvector3_t(a.x*k,a.y*k,a.z*k);
    }
  
  //& operator
  //vectorial product
  inline cvector3_t operator&(const vector3_t & a,const cvector3_t & b){
    return cvector3_t(a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x);
    }

class vector2D_t {
private:
public:
  double  x,y;

  vector2D_t(){
    x=0.;
    y=0.;
    }

  void init(double x0,double y0){
    x=x0;
    y=y0;
    }

  vector2D_t(double x0,double y0){
    init(x0,y0);
    }

  vector2D_t(point2D_t p, point2D_t q) {
    x=q.x-p.x;
    y=q.y-p.y;
    }

  vector2D_t(double alpha){
    x=cos(alpha);
    y=sin(alpha);
    }

  inline vector2D_t(double x1, double y1, double x2, double y2) {
    x=x2-x1;
    y=y2-y1;
    }

  inline vector2D_t& operator += (const vector2D_t& __x) {
    x +=  __x.x;
    y +=  __x.y;
    return *this;
    }

  inline vector2D_t& operator -= (const vector2D_t& __x) {
    x -=  __x.x;
    y -=  __x.y;
    return *this;
    }

  inline vector2D_t& operator *= (const double c) {
    x *=  c;
    y *=  c;
    return *this;
    }

  inline vector2D_t& operator /= (const double c) {
    x /=  c;
    y /=  c;
    return *this;
    }
  
  inline vector2D_t normal(){
    vector2D_t __r;
    double n=sqrt(x*x+y*y);
    __r.x=-y/n;
    __r.y= x/n;
    return __r;
    }

  };
  
  inline double square(const vector2D_t & v){
    double d2;
    d2=square(v.x)+square(v.y);
    return d2;
    }

  inline double hypot(const vector2D_t & v){
    double d;
    d=hypot(v.x,v.y);
    return d;
    }

  inline vector2D_t operator * (const double& __c, const vector2D_t& __x) {
    vector2D_t __r = __x;
    __r *= __c;
    return __r;
    }

  inline vector2D_t operator * (const vector2D_t& __x, const double __c) {
    vector2D_t __r = __x;
    __r *= __c;
    return __r;
    }

  inline vector2D_t operator / (const vector2D_t& __x, const double __c) {
    vector2D_t __r = __x;
    __r /= __c;
    return __r;
    }

  inline vector2D_t operator + (const vector2D_t& __x, const vector2D_t& __y) {
    vector2D_t __r = __x;
    __r += __y;
    return __r;
    }

  inline point2D_t operator += (const point2D_t& __x, const vector2D_t& __y) {
    point2D_t __r = __x;
    __r.x += __y.x;
    __r.y += __y.y;
    return __r;
    }

  inline vector2D_t operator - (const vector2D_t& __x, const vector2D_t& __y) {
    vector2D_t __r = __x;
    __r -= __y;
    return __r;
    }

/* *----------------------------------------------------------------------------
  2D vector translation */
  inline point2D_t operator + (const point2D_t& __x, const vector2D_t& __y) {
    point2D_t __r = __x;
//    __r += __y; /// HERE
    __r.x += __y.x;
    __r.y += __y.y;
    return __r;
    }

/* *----------------------------------------------------------------------------
  2D vector product */
  inline double operator * (const vector2D_t& __x, const vector2D_t& __y) {
    double __r;
    __r = __x.x*__y.y-__y.x*__x.y;
    return __r;
    }

/* *----------------------------------------------------------------------------
  2D scalar product */
  inline double operator % (const vector2D_t& __x, const vector2D_t& __y) {
    double __r;
    __r = __x.x*__y.x+__y.y*__x.y;
    return __r;
    }
  
/* *----------------------------------------------------------------------------
  2D angle */
  inline double operator^(const vector2D_t & a,const vector2D_t & b){
    double angle=atan2(a*b,a%b);
    return angle;
    }

class  cvector2D_t
{
private:
public:
  complex< double >  x,y;
  
  cvector2D_t(){
    x=0;
    y=0;
    }

  cvector2D_t(const complex<double> & x0,const complex<double> & y0){
    x=x0;
    y=y0;
    }

  };

  inline complex<double> operator * (const cvector2D_t& __x, const vector2D_t& __y) {
    complex< double >   __r;
    __r = __x.x*__y.y-__y.x*__x.y;
    return __r;
    }

  inline complex<double> operator * (const vector2D_t& __x, const cvector2D_t& __y) {
    complex< double >   __r;
    __r = __x.x*__y.y-__y.x*__x.y;
    return __r;
    }


  inline cvector2D_t operator * (const zmatrix2x2_t& __c, const cvector2D_t& __x) {
    cvector2D_t __r;
    __r.x = __c.c[0][0]*__x.x+__c.c[0][1]*__x.y;
    __r.y = __c.c[1][0]*__x.x+__c.c[1][1]*__x.y;
    return __r;
    }
  
extern matrix2x2_t math_rotation2D_init(double angle);
extern vector2D_t  math_rotation2D(const matrix2x2_t & R,const vector2D_t & u);

extern matrix3x3_t math_rotation3D_init(const vector3_t & axis, double angle);
extern matrix3x3_t math_rotation3D_init(char a, double angle);
extern vector3_t math_rotation3D(char operation, const matrix3x3_t & R,const vector3_t & u);

extern int math_rotation3Ddeg(char operation,const matrix3x3_t & R1,const matrix3x3_t & R2, double *x,double *y);
extern int math_rotation3Ddeg(char operation,const matrix3x3_t & R, double *x,double *y);

extern vector3_t math_polar2cartesian(double lon, double lat, double radius=1.0);
extern int math_cartesian2polar(const vector3_t & v,double *lon, double *lat, double *radius=NULL);

extern int math_orthobase(const point2D_t & p, const point2D_t & q, vector2D_t *u, vector2D_t *v, double *normp=NULL);

extern vector2D_t math_vector_coordinates01(const vector2D_t & r,const vector2D_t & u,const vector2D_t & v);
extern vector2D_t math_vector_coordinates02(const vector2D_t & r,const vector2D_t & u,const vector2D_t & v);

extern vector2D_t  math_vector_coordinates03(double a,const vector2D_t & u, double b,const vector2D_t & v);
extern cvector2D_t math_vector_coordinates03(complex<double> a,const vector2D_t & u, complex<double> b,const vector2D_t & v);

extern double hypot(double x,double y,double z);
extern double manhattan(double x,double y,double z=0.);

template<typename T> inline T square(const T& x){return x*x;}

template<typename T> T sign(T x){return (x < 0) ? -1 : (x > 0);}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> T operator%(const complex<T>& z1, const complex<T>& z2)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/// scalar product of 2 complex vectors
/// \note % has same precedence as * and /
{
  return real(z1) * real(z2) + imag(z1) * imag(z2);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> bool isnan(const complex<T>& z)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  bool result;
  
  result=isnan(real(z)) || isnan(imag(z));
  
  return result;
}

typedef struct {
  int x;
} monomial1D_t;

typedef struct {
  int x, y;
} monomial2D_t;

typedef struct {
  int x, y, z;
} monomial3D_t;

typedef struct {
  monomial2D_t **m;
  double **c;
  int order, n;
} polynomial_t;

typedef struct {
  monomial3D_t ***m;
  double ***c;
  int order, n;
} polynomial3D_t;



extern  int matrix_print(double *matrix, int nrow, int ncol);

extern int matrix_inverse(double *matrix, int neq, double* & inverse, bool preserve, int verbose);

extern int matrix_transpose(double *matrix, int nrow, int ncol, double* & transpose, int verbose);

extern int matrix_product(double *A, double *B, int nrowA, int ncolA, int nrowB, int ncolB, double* & product, int verbose);

extern int matrix_product(double *A, double *B, int ndim, double* & product, int verbose);


// template <typename P, typename Q> int matrix_operation(P *A, Q *B, int nrowA, int ncolA, Q* & product, int verbose)
extern int matrix_operation(double *A, double *B, int nrowA, int ncolA, double* & product, int verbose);
extern int matrix_operation(double *A, complex<double> *B, int nrowA, int ncolA, complex<double>* & product, int verbose);


extern int matrix_least_square(double *matrix, int nrow, int ncol, double* & leastsquare, int verbose);


#endif
