
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Yoann Le Bars      LEGOS, Toulouse, France
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada

***************************************************************************/

#ifndef __CEFMO
#define __CEFMO

#include "tools-define.h"
#include "tools-structures.h"
#include "tides.h"
#include <map>

class tide2D_t {
private :
public :
/*-----------------------------------------------------------------------------------------------------
  physical parameters*/
  double *h_z,*h_u;
  complex <double> *z,*u,*v;              /* wave drag rugosty, bottom friction coefficient, bottom rugosity*/
  complex <double> *LSA,*potential;       /* wave drag rugosty, bottom friction coefficient, bottom rugosity*/

  tide2D_t() {
    h_z=0;
    h_u=0;
    u=0;
    v=0;
    z=0;
    LSA=0;
    potential=0;
    }

  void destroy() {
    deletep(&h_z);
    deletep(&h_u);
    deletep(&u);
    deletep(&v);
    deletep(&z);
    deletep(&LSA);
    deletep(&potential);
    }
};

class tide3D_t {
private :
public :
/**----------------------------------------------------------------------------
  physical parameters*/
  double *h_z,*h_u;
  complex <double> **z,**u,**v/*,**w*/;    /* level displacement, layer velocities             */
  complex <double> *u_mean,*v_mean;    /* mean velocities                                  */
  complex <double> *LSA,*potential;    /* loading/self-attraction, astronomical potential  */
  double **rho;                        /* density                                          */
  int n_Unodes,n_Znodes,nlevels,nlayers; /*                                                  */

tide3D_t() {
    h_z=0;
    h_u=0;
    u=0;
    v=0;
//     w=0;
    u_mean=v_mean=0;
    z=0;
    LSA=0;
    potential=0;
    rho=0;
    n_Unodes=n_Znodes=nlayers=nlevels=-1;
    }
    
//   void allocate(int n1, int n2, int n3) {
//     if(n1<=0) return;
//     this->nnodes=n1;
//     u_mean=new complex <double> [nnodes];
//     v_mean=new complex <double> [nnodes];
//     if(n2>0) {
//       z=new complex <double> *[nnodes];
//       for(size_t n=0;n<nnodes;n++) {
//         z[n]=new complex <double> [n2];
//         }
//       this->nlevels=n2;
//       }
//     if(n3>0) {
//       u=new complex <double> *[nnodes];
//       v=new complex <double> *[nnodes];
//       for(size_t n=0;n<nnodes;n++) {
//         u[n]=new complex <double> [n3];
//         v[n]=new complex <double> [n3];
//         }
//       this->nlayers=n3;
//       }
//     }
    
  void allocate(int n1, int n2, int n3, int n4) {
    if(n1<=0) return;
    this->n_Unodes=n1;
    this->n_Znodes=n2;
    u_mean=new complex <double> [n_Unodes];
    v_mean=new complex <double> [n_Unodes];
    if((n2>0) && (n4>0)) {
      z  =new complex <double> *[n_Znodes];
      rho=new double *[n_Znodes];
      for(size_t n=0;n<n_Znodes;n++) {
        z[n]  =new complex <double> [n4];
        rho[n]=new double [n4];
        }
      this->nlevels=n4;
      }
    if(n3>0) {
      u=new complex <double> *[n_Unodes];
      v=new complex <double> *[n_Unodes];
      for(size_t n=0;n<n_Unodes;n++) {
        u[n]=new complex <double> [n3];
        v[n]=new complex <double> [n3];
        }
      this->nlayers=n3;
      }
    }
    
  void destroy() {
    deletep(&h_z);
    deletep(&h_u);
    deletep(&u_mean);
    deletep(&v_mean);
    deletep2D(&z,n_Znodes);
    deletep2D(&rho,n_Znodes);
    deletep2D(&u,n_Unodes);
    deletep2D(&v,n_Unodes);
//     deletep2D(&w);
    deletep(&LSA);
    deletep(&potential);
    }

  void duplicate() {//unused
    }
};



class actionZ_t {
private :
public :
  complex <double>  *prx, *pry;      /* pressure gradient                         */
  complex <double>  *tpx, *tpy;      /* potential                                 */
  complex <double>  *wdx, *wdy;      /* wave drag                                 */
  complex <double>  *Fx,  *Fy;       /* miscellaneous                             */

  actionZ_t() {
    prx=pry=0;
    tpx=tpy=0;
    wdx=wdy=0;
    Fx=Fy=0;
    }

  void destroy() {
    deletep(&prx);
    deletep(&pry);
    deletep(&tpx);
    deletep(&tpy);
    deletep(&wdx);
    deletep(&wdy);
    deletep(&Fx);
    deletep(&Fy);
    }

};

class cefmo_t {
private :
public :
/*-----------------------------------------------------------------------------------------------------
  physical parameters*/
  zmatrix2x2_t    *FrictionMatrix;
  hyperzmatrix_t   SpMatrix;
  matrix2x2_t     *DragMatrix;
// // // // // // //  double *R[5];
  int *clamped, nclamped;
  int paire;
  tide2D_t state;
  actionZ_t action;
  complex <double>  *NLRhs;
  tidal_wave dominants[2];
  int ndominants;

  cefmo_t() {
//    int k;
//    for(k=0;k<5;k++) R[k]=0;
    FrictionMatrix=0;
    DragMatrix=0;
    clamped=0;
    nclamped=-1;
    paire=-1;
    NLRhs=0;
    }

  void destroy() {
//    int k;
//    for(k=0;k<5;k++) deletep(&R[k]);
    deletep(&FrictionMatrix);
    deletep(&DragMatrix);
    deletep(&clamped);
    state.destroy();
    action.destroy();
    nclamped=-1;
    paire=-1;
    deletep(&NLRhs);
    }
};


typedef std::map<const char *, tide2D_t, cmp_str> atlas2D_t;
typedef std::map<const char *, tide3D_t, cmp_str> atlas3D_t;

extern void accelere(double rn0,double rn1,double *rn2);

extern  void momentum2d_DGP1xLGP2_SpInitialise(tidal_wave wave, cefmo_t model, parameter_t data,mesh_t mesh);
extern  void momentum2d_DGP1xLGP2_SpInitialise_iterate(tidal_wave wave,cefmo_t model,mesh_t mesh);

extern  complex<double>  *WE2D_DGP1xLGP2_SpInitialise(mesh_t mesh,tidal_wave wave,cefmo_t model,parameter_t data, ordering_t *ordering);
extern  void WE2D_DGP1xLGP2_SpInitialise_iterate(mesh_t mesh,tidal_wave wave,cefmo_t model,hyperzmatrix_t elevation);

extern  void WE2D_DGP1xLGP2_SpRHS(mesh_t mesh, cefmo_t model, parameter_t data, complex<double>  *rhs);
extern  void momentum_DGP1xLGP2_SpSolver(mesh_t mesh, cefmo_t model, parameter_t data);

extern void cefmo_frictionRHS(spectrum_t, spectrum_t, double,tide2D_t *, parameter_t, int, complex<double> **, complex<double> **);

extern zmatrix2x2_t *spectral_friction01(mesh_t &, complex <double> *, complex <double> *, int);
extern zmatrix2x2_t *spectral_friction02(mesh_t &, complex <double> *, complex <double> *, int);


extern void spectral_friction02(mesh_t, zmatrix2x2_t *, complex <double> *, complex <double> *, int);
extern void friction_coefficient(zmatrix2x2_t FrictionMatrix, double Cd, double r, double h, complex<double> *r1, complex<double> *r2, complex<double> *r3, complex<double> *r4);
extern void friction_coefficient(zmatrix2x2_t FrictionMatrix, double Cd, double r, double u0, double h, double ratio, complex<double> *r1, complex<double> *r2, complex<double> *r3, complex<double> *r4);


#endif