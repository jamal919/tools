
/**************************************************************************

  Double Kelvin Wave for canal test

  
Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  Yves Soufflet      LEGOS/CNRS, Toulouse, France

Date: 26/01/2011

***************************************************************************/

#define MAIN_SOURCE

#include <stdio.h>
#include <stdlib.h>
#include <cmath>
//#include "constants.h"
#include <complex>

#include "config.h"


#include "fe.h"

#include "geo.h"
#include "functions.h"
#include "tools-structures.h"
#include "tides.h"
#include "statistic.h"

#include "rutin.h"     /*  rutin.h contains common utility routines  */

using namespace std;

extern int get_openlimits(const char *,const char *, double **, double **, int *);


/**************************************************************************
Considering a channel of length L, and with l:


\eta(x,y)=K\left[\exp\left(jkx+fky\w^{-1}\right)+\exp\left(-jkx-fk(y-l)\w^{-1}\right)\right]

ele(x,y) = K exp(ikx+fky/omega) + K exp(-ikx-fk(y-l)/omega)

with a damping term q=lfk/(2.omega.lambda) and a closed boundary:

ele(x,y) = K exp((ik-q)x+fky/omega) + K exp((2L-x)(ik-q)-fk(y-l)/omega)

**************************************************************************/

// c-------------------------------------------------------------------------------------------------
// c
//       subroutine amphi(xlo,xla,a,p)
// c
// c   Calcul de la solution analytique pour deux ondes se deplacant
// c   en sens inverse dans un canal.
// c   Attention amplitudes divisees par deux et phases changes de signes
// c
// c-------------------------------------------------------------------------------------------------
// c
//
//       common/const1/rlamins,rlamaxs,rlomins,rlomaxs
//       common/const2/xt1,freq,xko,g,depth,f,theta,xk
//
//       nr=10
//       no=20
//
//
//       freq=1.4E-04
//       rlamins=50.0-0.899
//       rlamaxs=50.0+0.899
//       rlomins=+0.000
//       rlomaxs=+4.4108
//       depth=50.
//       r1=3.12E-05
//       ray=6.378E+06
//       pi=acos(-1.E+00)
//       pis180=pi/180.E+00
//       xk=3.1
//       g=9.81E+00
//       theta=.5*atan(r1/freq)
//       xko=sqrt(1.+(r1*r1)/(freq*freq))
//       xt1=pis180*ray
//       f=(pi/21600.)*sin(50.*pis180)
//
//       complex z
//
//       y=xt1*(xla-((rlamins+rlamaxs)/2.))
//       x=xt1*(xlo-rlomins)
//       xx=x*freq*sqrt(xko/(g*depth))
//       yy=y*sqrt(1./(g*depth*xko))*f
//       a1=yy*sin(theta)+xx*cos(theta)
//       a2=yy*cos(theta)+xx*sin(theta)
//       v=sin(a1)*(exp(a2)-1.*exp(-a2))
//       u=cos(a1)*(exp(a2)+1.*exp(-a2))
//       u=xk*u
//       v=(-1.)*xk*v
//
//       z=cmplx(u,v)
//       a=ro(z)
//       p=beta(z)
//
//       return
//       end

//Global variables
double g=9.81; // Earth's gravity
#define pi M_PI
complex<double> I(0.,1.);

class kelvin_wave
{
private:
public:
    double lati,L,l,h,T,omega,K;

    kelvin_wave() {
        printf("in constructor, setting values \n");
        double lati=0.; // latitude
        double L= 484992.791; //Length of basin 750km
        double l= 197708.8538; //width of basin 250 km  L=750E3;
        double h=40.; // depth that gives the most realistic response

        // M2 tidal wave
        double T=12.75*3600; // period
        double omega=2*pi/T; // pulsation
        double K=pi; // arbitrary tidal amplitude. pi is good for isolines because it macthes max phase values
    }

    void print(void)
    {
        printf("lati=%lf \n",lati);
        //std::cout<<"lati="<<lati<<"\n";
        printf("L=%lf \n",L);
        printf("l=%lf \n",l);
        printf("T=%lf \n",T);
        printf("omega=%lf \n",omega);
        printf("K=%lf \n",K);
        //std::cout<<"K="<<K<<"\n";
    }

    int double_wave_open(double *lon,double *lat,int count, double *Ha,double *Hg)
    {
        //sale boulot en attendant fix
        double lati=45.; // latitude
        double L= 484992.791; //Length of basin 750km
        double l= 197708.8538; //width of basin 250 km  L=750E3;
        double h=40.; // depth that gives the most realistic response

        // M2 tidal wave
        double T=12.75*3600; // period
        double omega=2*pi/T; // pulsation
        double K=pi; // arbitrary tidal amplitude. pi is good for isolines because it macthes max phase values
       
       /*-----------------------------------------------------------------------------
         work in cartesian coordinates */
       double radius=6300.0,angle;
        double *x,*y;
        geo_t projection;
        int status=1;
        printf("g=%f \n",g);
        printf("h=%f \n",h);
        x=(double *) malloc(count*sizeof(double));
        y=(double *) malloc(count*sizeof(double));
        printf("applying the boundary conditions to cartesian coordinates\n");
        double f=two_Omega*sin(lati); // Coriolis parameter
        printf("Coriolis parameter= %f \n",f);
        double c=sqrt(g*h); // wave velocity
        printf("wave velocity=%f \n",c);
        double k=omega/c; // wave number
        printf("wave number=%f \n",k);
        double lambda=c*T; // wave length (lambda in descritption)
        double q=l*f*k/(omega*2*lambda);//dumping term
        complex<double> eta;
        projection=geo_mercator_init(0.0,0.0,radius);
        for (int n=0;n<count;n++) {

            status=geo_mercator_directe(projection,lon[n],lat[n],&x[n],&y[n]);
            printf("in double wave, values lon lat are: %f %f\n",lon[n],lat[n]);
            printf("in double wave, values x,y are: %f %f\n",x[n],y[n]);

            printf("imaginary and real parts are: %f %f\n", k*x[n],f*k*y[n]/omega);
            printf("imaginary and real parts are: %f %f\n", -1.0*k*x[n],-1.0*f*k*(y[n]-l)/omega);
            eta=K*exp(I*k*x[n]+f*k*y[n]/omega)+K*exp(-1.0*I*k*x[n]-f*k*(y[n]-l)/omega);
            Ha[n]=abs(eta);
            Hg[n]=arg(eta);
            printf("in double wave, values are: %f %f\n",Ha[n],Hg[n]);
        }
        return status=0;
    }

    complex<double> double_wave_semiclosed(double *lon,double *lat,int count, double *Ha,double *Hg)
    {//sale boulot en attendant fix
        double lati=45.; // latitude
        double L= 484992.791; //Length of basin 750km
        double l= 197708.8538; //width of basin 250 km  L=750E3;
        double h=40.; // depth that gives the most realistic response

        // M2 tidal wave
        double T=12.75*3600; // period
        double omega=2*pi/T; // pulsation
        double K=pi; // arbitrary tidal amplitude. pi is good for isolines because it macthes max phase values
       
       /*-----------------------------------------------------------------------------
         work in cartesian coordinates */
       double radius=6300.0,angle;
        double *x,*y;
        geo_t projection;
        int status=1;
        printf("g=%f \n",g);
        printf("h=%f \n",h);
        x=(double *) malloc(count*sizeof(double));
        y=(double *) malloc(count*sizeof(double));
        printf("applying the boundary conditions to cartesian coordinates\n");
        double f=two_Omega*sin(lati); // Coriolis parameter
        printf("Coriolis parameter= %f \n",f);
        double c=sqrt(g*h); // wave velocity
        printf("wave velocity=%f \n",c);
        double k=omega/c; // wave number
        printf("wave number=%f \n",k);
        double lambda=c*T; // wave length (lambda in descritption)
        double q=l*f*k/(omega*2*lambda);//dumping term
        complex<double> eta;
        projection=geo_mercator_init(0.0,0.0,radius);
        for (int n=0;n<count;n++) {

            status=geo_mercator_directe(projection,lon[n],lat[n],&x[n],&y[n]);
            printf("in double wave, values lon lat are: %f %f\n",lon[n],lat[n]);
            printf("in double wave, values x,y are: %f %f\n",x[n],y[n]);

            printf("imaginary and real parts are: %f %f\n", k*x[n],f*k*y[n]/omega);
            printf("imaginary and real parts are: %f %f\n", -1.0*k*x[n],-1.0*f*k*(y[n]-l)/omega);
            eta=K*exp(I*k*x[n]+f*k*y[n]/omega)+K*exp(-1.0*I*k*x[n]-f*k*(y[n]-l)/omega);
            Ha[n]=abs(eta);
            Hg[n]=arg(eta);
            printf("in double wave, values are: %f %f\n",Ha[n],Hg[n]);
        }
        return status=0;
    }

    int provost_vincent(double *lon, double *lat, int count, double *Ha,double *Hg, double *Ua,double *Ug, double *Va, double *Vg)
      {
        //sale boulot en attendant fix
        
/*-----------------------------------------------------------------------------
        Le Provost, C., and P. Vincent, Some tests of precision for a finite
        element model of ocean tides, J. Comput. Phys., 65, 273-291, 1986. */

        double lati=45.;        // latitude
        double L= 484992.791;   // Length of basin 750km
        double l= 197708.8538;  // width of basin 250 km  L=750E3;
        double h=50.;           // depth
        
        // M2 tidal wave
        //double T=12.75*3600;  // period
        double omega=1.4e-04;   // pulsation
        double K=3.0;           // arbitrary tidal amplitude.
        double r=3.12e-05;
        double k0=sqrt(1+r*r/(omega*omega));
        double theta=0.5*atan(r/omega);
        
/*-----------------------------------------------------------------------------
        work in cartesian coordinates */
        double radius=6300.0,angle;
        double *x,*y;
        double X,Y;
        geo_t projection;
        int status=1;
        printf("g=%f \n",g);
        printf("h=%f \n",h);
        x= new double[count];
        y= new double[count];
        double f=two_Omega*sin(lati); // Coriolis parameter
        complex<double> eta;
        complex<double> mu;
/*-----------------------------------------------------------------------------
        compute cartesian coordinates */
        projection=geo_mercator_init(0.0,0.0,radius);
        for (int n=0;n<count;n++) {
          status=geo_mercator_directe(projection,lon[n],lat[n],&x[n],&y[n]);
          }
        range_t<double> rangeX=poc_minmax(x,count,1e10);
        range_t<double> rangeY=poc_minmax(y,count,1e10);
        for (int n=0;n<count;n++) {
          x[n]-=rangeX.min;
          y[n]-=rangeY.min;
          }
        rangeX=poc_minmax(x,count,1e10);
        rangeY=poc_minmax(y,count,1e10);
/*-----------------------------------------------------------------------------
        compute analytic solution */
        for (int n=0;n<count;n++) {
//          printf("in double wave, values lon lat are: %f %f\n",lon[n],lat[n]);
          printf("x,y : %lf %lf\n",x[n],y[n]);
          
          X=x[n]*omega*sqrt(k0/(g*h));
          Y=y[n]*f/sqrt(g*h*k0);
            
          eta=K*exp(X*cos(theta)+Y*sin(theta))*exp(I*(X*sin(theta)+Y*cos(theta)));

//  eta=K*exp(I*k*x[n]+f*k*y[n]/omega)+K*exp(-1.0*I*k*x[n]-f*k*(y[n]-l)/omega);

          Ha[n]=abs(eta);
          Hg[n]=arg(eta);
          printf("in double wave, values are: %f %f\n",Ha[n],Hg[n]);
          mu=-1.0*sqrt(k0/(g*h))*K*exp(I*theta)*eta;
          Ua[n]=abs(mu);
          Ug[n]=arg(mu);
          Va[n]=0.0;
          Vg[n]=0.0;
          }
        return status=0;
      }

    int provost_vincent(discretisation_t z_descriptor, discretisation_t u_descriptor, complex<double> *z, complex<double> *u, complex<double> *v)
      {
        //sale boulot en attendant fix
        
/*-----------------------------------------------------------------------------
        Le Provost, C., and P. Vincent, Some tests of precision for a finite
        element model of ocean tides, J. Comput. Phys., 65, 273-291, 1986. */

//       freq=1.4E-04

//       rlamins=50.0-0.899
//       rlamaxs=50.0+0.899
//       rlomins=+0.000
//       rlomaxs=+4.4108

//       depth=50.
//       r1=3.12E-05
//       ray=6.378E+06

//       pis180=pi/180.E+00

//       A=3.1

//       g=9.81E+00
//       theta=.5*atan(r1/freq)
//       xko=sqrt(1.+(r1*r1)/(freq*freq))
//       xt1=pis180*ray
//       f=(pi/21600.)*sin(50.*pis180)

//       y=xt1*(xla-((rlamins+rlamaxs)/2.))
//       x=xt1*(xlo-rlomins)

//       X=x*freq*sqrt(xko/(g*depth))
//       Y=y*sqrt(1./(g*depth*xko))*f

//       a1=Y*sin(theta)+X*cos(theta)
//       a2=Y*cos(theta)+X*sin(theta)

//       u=cos(a1)*(exp(a2)+1.*exp(-a2))
//       v=sin(a1)*(exp(a2)-1.*exp(-a2))

//       eta1=exp( a2) * (cos(a1),  sin(a1)) = exp( a2) * exp(+J*a1)
//       eta2=exp(-a2) * (cos(a1), -sin(a1)) = exp(-a2) * exp(-J*a1)

//       u=+A*u
//       v=-A*v
//
//       z=cmplx(u,v)
//       a=ro(z)
//       p=beta(z)
        double lati=50.;        // latitude
        double L= 484992.791;   // Length of basin 750km
        double l= 197708.8538;  // width of basin 250 km  L=750E3;
        double h=50.;           // depth
        
        // M2 tidal wave
        //double T=12.75*3600;  // period
        double omega=1.4e-04;   // pulsation
        double T=2*M_PI/omega;     // period
        double K=3.0;           // arbitrary tidal amplitude.
        double r=3.12e-05;
/*-----------------------------------------------------------------------------
        Kelvin wave decay scale (adimensional) */
        double decay=r/omega;
        double k0=sqrt(1+decay*decay);
        double theta=0.5*atan(decay);
        double c=sqrt(g*h);
        double lambda=c*T, amplitude;
        
/*-----------------------------------------------------------------------------
        work in cartesian coordinates */
        double radius=6300.0,angle;
        double *x,*y;
        double *lon, *lat;
        double X,Y;
        double *Ha, *Hg, *Ua, *Ug, *Va, *Vg;

        geo_t projection;
        int status=1;
        
        range_t<double> rangeX, rangeY;
        
        printf("g=%f \n",g);
        printf("h=%f \n",h);

        double f=two_Omega*sin(lati*d2r); // Coriolis parameter
        complex<double> eta;
        complex<double> mu;

        projection=geo_mercator_init(0.0,0.0,radius);

/*-----------------------------------------------------------------------------
        compute cartesian coordinates */
        x= new double[z_descriptor.nnodes];
        y= new double[z_descriptor.nnodes];
        for (int n=0;n<z_descriptor.nnodes;n++) {
          status=geo_mercator_directe(projection,z_descriptor.nodes[n].lon,z_descriptor.nodes[n].lat,&x[n],&y[n]);
          }

/*-----------------------------------------------------------------------------
        normalize cartesian coordinates (to fit publication parameters) */
        for (int n=0;n<z_descriptor.nnodes;n++) {
//           x[n]*=1000./lambda;
//           y[n]*=1000./lambda;
          x[n]*=1000.;
          y[n]*=1000.;
          }
        rangeX.init(x,z_descriptor.nnodes);
        rangeY.init(y,z_descriptor.nnodes);
        for (int n=0;n<z_descriptor.nnodes;n++) {
          x[n]-=rangeX.min;
          y[n]-=rangeY.min;
          }
        rangeX.init(x,z_descriptor.nnodes);
        rangeY.init(y,z_descriptor.nnodes);
        for (int n=0;n<z_descriptor.nnodes;n++) {
          x[n]*=500000./rangeX.max;
          y[n]*=200000./rangeY.max;
          y[n]-=100000.;
          }
        rangeX.init(x,z_descriptor.nnodes);
        rangeY.init(y,z_descriptor.nnodes);

/*-----------------------------------------------------------------------------
        compute analytic solution */
        for (int n=0;n<z_descriptor.nnodes;n++) {
//          printf("in double wave, values lon lat are: %f %f\n",lon[n],lat[n]);
          printf("x,y : %lf %lf\n",x[n],y[n]);

/*-----------------------------------------------------------------------------
          X=x*(2*pi)/lamda(wave)*sqrt(k0) */
          X=x[n]*omega*sqrt(k0)/c;
/*-----------------------------------------------------------------------------
          Y=y*(2*pi)/lamda(inertial)/sqrt(k0) */
          Y=y[n]*f    /sqrt(k0)/c;

          eta=0.0;

/*-----------------------------------------------------------------------------
          eastward incoming wave (from western boundary) */
          amplitude=K*exp(-(Y*cos(theta)+X*sin(theta)));
          eta+=amplitude*exp(-I*(Y*sin(theta)+X*cos(theta)));

          X=(-x[n])*omega*sqrt(k0)/c;
          Y=(-y[n])*f    /sqrt(k0)/c;

/*-----------------------------------------------------------------------------
          westward incoming wave (from eastern boundary) */
          amplitude=K*exp(-(Y*cos(theta)+X*sin(theta)));
          eta+=amplitude*exp(-I*(Y*sin(theta)+X*cos(theta)));

          z[n]=eta;
          }
        delete[] x;
        delete[] y;

/*-----------------------------------------------------------------------------
        compute cartesian coordinates */
        x= new double[u_descriptor.nnodes];
        y= new double[u_descriptor.nnodes];
        for (int n=0;n<u_descriptor.nnodes;n++) {
          status=geo_mercator_directe(projection,u_descriptor.nodes[n].lon,u_descriptor.nodes[n].lat,&x[n],&y[n]);
          }
        rangeX.init(x,u_descriptor.nnodes);
        rangeY.init(y,u_descriptor.nnodes);
        for (int n=0;n<u_descriptor.nnodes;n++) {
          x[n]-=rangeX.min;
          y[n]-=rangeY.min;
          }
        rangeX.init(x,u_descriptor.nnodes);
        rangeY.init(y,u_descriptor.nnodes);

/*-----------------------------------------------------------------------------
        compute analytic solution */
        for (int n=0;n<u_descriptor.nnodes;n++) {
//          printf("in double wave, values lon lat are: %f %f\n",lon[n],lat[n]);
          printf("x,y : %lf %lf\n",x[n],y[n]);
          X=x[n]*omega*sqrt(k0/(g*h));
          Y=y[n]*f/sqrt(g*h*k0);
            
          eta=0.0;
          eta+=K*exp(Y*cos(theta)+X*sin(theta))*exp(I*(Y*sin(theta)+X*cos(theta)));
          X=-X;
          Y=-Y;
          eta+=K*exp(Y*cos(theta)+X*sin(theta))*exp(I*(Y*sin(theta)+X*cos(theta)));
          mu=-1.0*sqrt(k0/(g*h))*K*exp(I*theta)*eta;
          u[n]=mu;
          v[n]=0.0;
          }

        delete[] x;
        delete[] y;

        return status=0;
    }
};


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int generate_obc(kelvin_wave kw, char *belfile,char *meshfile, char *output, char **wave, int nwave)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int *code;
  int i,k,l,m,n,n1,n2,count,nndes,status;
  int nopen=0;
  double **Ha,**Hg,**Ua,**Ug,**Va,**Vg,mask;
  double a[3],G[3];
  double *lon,*lat;
  char filename[1024];
  FILE *out;

  status=get_openlimits(belfile,meshfile,&lon,&lat,&count);

  Ha=new double*[nwave];
  Hg=new double*[nwave];
  Ua=new double*[nwave];
  Ug=new double*[nwave];
  Va=new double*[nwave];
  Vg=new double*[nwave];
  
  
  for (k=0;k<nwave;k++) {
    Ha[k]=new double[count];
    Hg[k]=new double[count];
    Ua[k]=new double[count];
    Ug[k]=new double[count];
    Va[k]=new double[count];
    Vg[k]=new double[count];
    }
 
  for (k=0;k<nwave;k++) {
/* *-----------------------------------------------------------------------------
    apply bc from function defined above */
    printf("applying wave bc %s wave \n",wave[k]);
    status=kw.double_wave_open(lon,lat,count, Ha[k],Hg[k]);
//    status=kw.provost_vincent(lon,lat,count, Ha[k],Hg[k],Ua[k],Ug[k],Va[k],Vg[k]);
//    status=tide_atlasSG2positions(atlas_file,"Ha","Hg",lon,lat,count,Ha[k],Hg[k],mask);
    if(status!=0) goto error;
/* *-----------------------------------------------------------------------------
    tidal heights and currents at open limits*/
    sprintf(filename,"%s.obc",wave[k]);
    out=fopen(filename,"w");
    fprintf(out,"%s\n",wave[k]);
    for(i=0; i<count; i++) {
/*---------------------------------------------------------------------
      elevation*/
      if(Hg[k][i]<0.0) Hg[k][i]+=360.0;
      if(Ug[k][i]<0.0) Ug[k][i]+=360.0;
      if(Vg[k][i]<0.0) Vg[k][i]+=360.0;
      a[0]=Ha[k][i];
      G[0]=Hg[k][i];
      //if((status_u!=0) || (status_v!=0)) {
        fprintf(out,"%8f %8f %f %f\n",lon[i],lat[i],a[0],G[0]);
        //}
      //else {
        //fprintf(out,"%8f %8f %f %f %f %f %f %f\n",lon[i],lat[i],a[0],G[0],Ua[k][i],Ug[k][i],Va[k][i],Vg[k][i]);
        //}
      }
    fclose(out);
    }

  if (output==0) {
    output=strdup("tides.obc");
    }

  out=fopen(output,"w");
  fprintf(out,"%d %d %s\n",count,nwave,"M");
  for (k=0;k<nwave;k++) {
/* *-----------------------------------------------------------------------------
    tidal heights at open limits*/
    fprintf(out,"%s\n",wave[k]);
    for(i=0; i<count; i++) {
/*---------------------------------------------------------------------
     elevation*/
      a[0]=Ha[k][i];
      G[0]=Hg[k][i];
      fprintf(out,"%8f %8f %f %f\n",lon[i],lat[i],a[0],G[0]);
      //fprintf(out,"%8f %8f %f %f %f %f %f %f\n",lon[i],lat[i],a[0],G[0],Ua[k][i],Ug[k][i],Va[k][i],Vg[k][i]);
      }
    }
  fclose(out);
  return(0);
  
error:
  printf("error detected, quit ... \n");
  return(-1);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int generate_solution2D(kelvin_wave kw, discretisation_t z_descriptor,discretisation_t  u_descriptor, complex<double> *z, complex<double> *u, complex<double> *v)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  
  for (n=0;n<z_descriptor.nnodes; n++) {
    }

}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int generate_solution2D(kelvin_wave kw, mesh_t& mesh, char *output, char **wave, int nwave, int paire)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int n,status;
  discretisation_t z_descriptor, u_descriptor;
  complex<double> *z,*u,*v, mask;
  
  paire=DNP1xLGP2;
  mesh.initialize_descriptors(paire, 0);
  
  z_descriptor=mesh.LGP2descriptor;
  u_descriptor=mesh.DNP1descriptor;
  
  z=new complex<double> [z_descriptor.nnodes];
  u=new complex<double> [u_descriptor.nnodes];
  v=new complex<double> [u_descriptor.nnodes];
  
  status=kw.provost_vincent( z_descriptor,  u_descriptor, z, u, v);

  status=archiving_UGdummy2D((const char*) output, mesh, "a_eta_LGP2", "G_eta_LGP2", "m",   z, mask, LGP2);
  status=archiving_UGdummy2D((const char*) output, mesh, "a_u_DNP1",   "G_u_DNP1",   "m/s", u, mask, DNP1);
  status=archiving_UGdummy2D((const char*) output, mesh, "a_v_DNP1",   "G_v_DNP1",   "m/s", v, mask, DNP1);

  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

int main(int argc, char *argv[])

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int count,nndes,status,fmt;
  int status_u,status_v;
  int i,j,k,l,m,n;
  int nitems;
  char *keyword;
  char *meshfile=NULL,*belfile=NULL,*output=NULL,*path=NULL,*atlas_convention=NULL;
  char *atlas_directory=NULL,*atlas_file,*format=0;
  char *wave[1024];
  fcomplex *z,*u,*v;
  mesh_t mesh;
  int nwave=0;
  spectrum_t spectrum;
  string s,executable,echofile;
  int pos;
  FILE *echo;
  kelvin_wave kw;
  fct_echo( argc, argv);
  
  n=1;
  while (n < argc) {
    keyword=strdup(argv[n]);
    switch (keyword[0]) {
      case '-':
        switch (keyword[1]) {
        case 'm' :
          meshfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'b' :
          belfile= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        naming convention for tidal atlases*/
        case 'c' :
          atlas_convention= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'f' :
          format= strdup(argv[n+1]);
          n++;
          n++;
          break;

        case 'o' :
          output= strdup(argv[n+1]);
          n++;
          n++;
          break;

/* *----------------------------------------------------------------------
        path for tidal atlases*/
        case 'p' :
          path= strdup(argv[n+1]);
          n++;
          n++;
          break;

        default:
          __OUT_BASE_LINE__("unknown option %s\n",keyword);
          exit(-1);
        }
        break;

      default:
        wave[nwave]= strdup(argv[n]);
        printf("input wave=%s\n",wave[nwave]);
        spectrum.n=nwave+1;
        nwave++;
        n++;
        break;
      }
      free(keyword);
    }

  if(path!=NULL) atlas_directory=strdup(path);

  if(format==0) format=strdup("T-UGOm");

  spectrum.waves=new tidal_wave[spectrum.n];
  for (k=0;k<nwave;k++) {
    strcpy(spectrum.waves[k].name,wave[k]);
    }

  if(path ==NULL) path=strdup(".");

  if(meshfile != NULL) {
    status=fe_readmesh(meshfile,MESH_FILE_FORMAT_TRIGRID,&mesh);
    if(status !=0) goto error;
    status=fe_list(&mesh);
    if(status !=0) goto error;
    status= fe_edgetable(&mesh,0,0);
    if(status !=0) goto error;
    }
  else {
    printf("no mesh file specified; abort...\n");
    goto error;
    }

//  generate_obc(kw,belfile, meshfile,  output,  wave,  nwave);

  status=generate_solution2D(kw, mesh, output, wave, nwave, DNP1xLGP2);

  printf("free memory \n");

end: __OUT_BASE_LINE__("end of %s ... \n",argv[0]);
  exit(0);
  
error:
 __OUT_BASE_LINE__("error detected, quit ... \n");
  exit(-1);
}



/*int main(){
  char meshfile[1024]="ploumploumploum";
  mesh_t mesh;
  
 kelvin_wave kv ;//= new kelvin_wave();
 printf("In main function of kelvin.cpp");
 complex<double> n;
 n=kv.double_wave_open_bc(0.0,0.0);
 printf("result is: %lf \n", abs(n));
 return(0);
}
 */
   
 
