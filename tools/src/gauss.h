
#define GAUSS_LGP0 0
#define GAUSS_LGP1 1
#define GAUSS_LGP2 2

#define GAUSS_NCP1 3

#define NGAUSS 4

class gauss_t {
private :
  int translation[1000];
  int discretisation[NGAUSS];
public :
  double *x, *y;
  double *w;
  double **base[NGAUSS],**base_x[NGAUSS],**base_y[NGAUSS];
  int nnodes,nnpe[NGAUSS];

  gauss_t() {
    int l;
    x=y=0;
    w=0;
    for (l=0;l<3;l++) {
      base[l]=base_x[l]=base_y[l]=0;
      nnpe[l]=-1;
      }
    nnodes=-1;
    translation[LGP0]=GAUSS_LGP0;
    translation[LGP1]=GAUSS_LGP1;
    translation[DGP1]=GAUSS_LGP1;
    translation[LGP2]=GAUSS_LGP2;
    translation[DGP2]=GAUSS_LGP2;
    translation[NCP1]=GAUSS_NCP1;
    translation[DNP1]=GAUSS_NCP1;
    discretisation[GAUSS_LGP0]=LGP0;
    discretisation[GAUSS_LGP1]=LGP1;
    discretisation[GAUSS_LGP2]=LGP2;
    discretisation[GAUSS_NCP1]=NCP1;
    }
  void destroy(){
    int k,l;
    deletep(&x);
    deletep(&y);
    deletep(&w);
    for (l=0;l<NGAUSS;l++) {
      deletep2D(&base[l],nnodes);
      deletep2D(&base_x[l],nnodes);
      deletep2D(&base_y[l],nnodes);
      }
    for (l=0;l<NGAUSS;l++) {
      nnpe[l]=-1;
      }
    nnodes=-1;
    }

  void init(int type){
    int k,l,n,status;
    mesh_t mesh;
    switch(type) {
      case IPG_7:
        n=7;
        break;
      }
    nnodes=n;
    x=new double[n];
    y=new double[n];
    w=new double[n];
    gauss_init(nnodes,x,y,w);
    nnodes=n;
    nnpe[GAUSS_LGP0]=1;
    nnpe[GAUSS_LGP1]=3;
    nnpe[GAUSS_LGP2]=6;
    nnpe[GAUSS_NCP1]=3;
    for (l=0;l<NGAUSS;l++) {
      base[l]  =new double*[nnodes];
      base_x[l]=new double*[nnodes];
      base_y[l]=new double*[nnodes];
      for (k=0;k<nnodes;k++) {
        base[l][k]  =new double[nnpe[l]];
        base_x[l][k]=new double[nnpe[l]];
        base_y[l][k]=new double[nnpe[l]];
/**----------------------------------------------------------------------------
        return interpolation function coefficients at integration points */
//         status=fe_LGPbase (mesh, discretisation[l], x[k], y[k], base[l][k]);
//         status=fe_LGPprime(mesh, discretisation[l], x[k], y[k], base_x[l][k], base_y[l][k]);
        }
      }
    }
    
  void coordinates(mesh_t mesh, int m, double *t, double *p){
    int k,i,n,discretisation,status;
    double tref;
        
//     n=mesh.triangles[m].vertex[0];
//     tref=mesh.vertices[n].lon;
//     for (k=0;k<nnodes;k++) {
//       t[k]=0;
//       p[k]=0;
//       for(i=0;i<3;i++) {
//         n=mesh.triangles[m].vertex[i];
//         t[k]+=base[GAUSS_LGP1][k][i]*geo_recale(mesh.vertices[n].lon,tref,180.0);
//         p[k]+=base[GAUSS_LGP1][k][i]*mesh.vertices[n].lat;
//         }
//       }
    }
    
template <typename T> T interpolate_template(T *p, int k, int discretisation){
    int i,id;
    T sum=0.0;
    id=translation[discretisation];
    for(i=0;i<nnpe[id];i++) {
      sum+=(T) base[id][k][i]*p[i];
      }
    return(sum);
    }
    
  double interpolate(float *p, int k, int discretisation){
    return(interpolate_template(p,k,discretisation));
    }
    
  double interpolate(double *p, int k, int discretisation){
    return(interpolate_template(p,k,discretisation));
    }
    
  complex<double> interpolate(complex<double> *p, int k, int discretisation){
    return(interpolate_template(p,k,discretisation));
    }
    
  complex<float> interpolate(complex<float> *p, int k, int discretisation){
    return(interpolate_template(p,k,discretisation));
    }
    
template <typename T> T interpolate_x_template(T *p, int k, int discretisation){
    int i,id;
    T sum=0.0;
    id=translation[discretisation];
    for(i=0;i<nnpe[id];i++) {
      sum+=base_x[id][k][i]*p[i];
      }
    return(sum);
    }
    
  double interpolate_x(float *p, int k, int discretisation){
    return(interpolate_x_template(p,k,discretisation));
    }
    
  double interpolate_x(double *p, int k, int discretisation){
    return(interpolate_x_template(p,k,discretisation));
    }
    
  complex<double> interpolate_x(complex<double> *p, int k, int discretisation){
    return(interpolate_x_template(p,k,discretisation));
    }
    
template <typename T> T interpolate_y_template(T *p, int k, int discretisation){
    int i,id;
    T sum=0.0;
    id=translation[discretisation];
    for(i=0;i<nnpe[id];i++) {
      sum+=base_y[id][k][i]*p[i];
      }
    return(sum);
    }
    
  double interpolate_y(float *p, int k, int discretisation){
    return(interpolate_y_template(p,k,discretisation));
    }
    
  double interpolate_y(double *p, int k, int discretisation){
    return(interpolate_y_template(p,k,discretisation));
    }
    
  complex<double> interpolate_y(complex<double> *p, int k, int discretisation){
    return(interpolate_y_template(p,k,discretisation));
    }
    
template <typename T> T integrate_template(T *p, int k, int discretisation){
    int i,id;
    T sum=0.0;
    id=translation[discretisation];
    for(i=0;i<nnpe[id];i++) {
      sum+=base[id][k][i]*p[i];
      }
    return(sum);
    }
    
};


extern float fe_QIntegraleLGP2xIPG7_2D(int i, float *q, gauss_t gauss);
  
extern double fe_QIntegraleLGP2xIPG7_2D(int i, double *q, gauss_t gauss);
  
extern complex<double> fe_QIntegraleLGP2xIPG7_2D(int i, complex<double> *q, gauss_t gauss);

