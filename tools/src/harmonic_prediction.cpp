
/*******************************************************************************

  T-UGO tools, 2006-2011

  Unstructured Ocean Grid initiative

*******************************************************************************/
/** \file

\author  Florent Lyard      LEGOS/CNRS, Toulouse, France. florent.lyard@legos.obs-mip.fr
\author  Laurent Roblou     LEGOS/CNRS, Toulouse, France
\author  Damien Allain      LEGOS/CNRS, Toulouse, France
\author  Yves Soufflet      LEGOS/CNRS, Toulouse, France
\author  Clément Mayet      LEGOS, Toulouse, France (PhD)
\author  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
\author  Frédéric Dupont    Université de Laval à Québec, Canada

\brief tidal prediction definitions
*/
/*----------------------------------------------------------------------------*/

#include "tides.h"
#include "poc-time.h"

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_predictionTS(double *serie,const double *time, int n,const spectrum_t & s, float *a, float *G, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Perform a tidal prediction for a single time serie
/*----------------------------------------------------------------------------*/
{
  int k,m;
  tidal_wave wave;
  double V,V0;
  double f,Vu;
  date_t start;
  astro_angles_t astro_angles;

  for(m = 0; m < n; m++) {
    serie[m] = 0.0;
    }

  for(k = 0; k < s.n; k++) {
    wave = s.waves[k];
    for(m = 0; m < n; m++) {
      start = poctime_getdatecnes(time[m], 'd');
      init_argument(&astro_angles,start);
      V0 = greenwhich_argument(astro_angles,wave);
      if(nodal == 1) {
        Vu = nodal_phase(astro_angles,wave);
        f = nodal_factor(astro_angles,wave.formula);
        }
      else {
        Vu = 0.;
        f = 1.;
        }
      V = V0 + Vu - G[k] * d2r;
      serie[m] += f * a[k] * cos(V);
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int harmonic_predictionTS_template_unsafe(T *serie, double *time, int n, spectrum_t s, fcomplex *constants, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Perform a tidal prediction for a single time serie
/*----------------------------------------------------------------------------*/
{
  int k,m;
  tidal_wave wave;
  double f,Vu;
  date_t start;
  extern date_t cnesdate(double t, char fmt);
  astro_angles_t astro_angles;

  for(m=0; m<n; m++) {
    serie[m]=0.0;
    }

  start=poctime_getdatecnes(time[0], 'd');
  init_argument(&astro_angles,start);

  for(k=0; k < s.n; k++) {
    wave=s.waves[k];
    double omega=wave.omega*dph2rpd;
    double V0=greenwhich_argument(astro_angles,wave);
    if(nodal==1) {
      Vu=nodal_phase(astro_angles,wave);
      f=nodal_factor(astro_angles,wave.formula);
      }
    else {
      Vu=0.;
      f=1.;
      }
#pragma omp parallel for
    for(m=0; m<n; m++) {
      double V=omega*(time[m]-time[0])+V0+Vu;
//       astro_angles_t astro_angles;
//       double V0=greenwhich_argument(astro_angles,wave);
//       if(nodal==1) {
//         Vu=nodal_phase(astro_angles,wave);
//         f=nodal_factor(astro_angles,wave.formula);
//         }
//       else {
//         Vu=0.;
//         f=1.;
//         }
//       double V=V0+Vu;
       
      T cs=cos(V);
      T sn=sin(V);
      serie[m]+=f*(real(constants[k])*cs-imag(constants[k])*sn);
      }
    }
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template <typename T> int harmonic_predictionTS_template(T *serie, double *time, int n, spectrum_t s, fcomplex *constants, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
/// Perform a tidal prediction for a single time serie
/*----------------------------------------------------------------------------*/
{
  int k,m;
  tidal_wave wave;
  extern date_t cnesdate(double t, char fmt);
  astro_angles_t *astro_angles;

  for(m=0; m<n; m++) {
    serie[m]=0.0;
    }

  astro_angles =new astro_angles_t[n];
  
#pragma omp parallel for
  for(m=0; m<n; m++) {
    date_t start=poctime_getdatecnes(time[m], 'd');
    init_argument(&astro_angles[m],start);
    }
    
  for(k=0; k < s.n; k++) {
    wave=s.waves[k];
#pragma omp parallel for
    for(m=0; m<n; m++) {
      double f,Vu;
      double V0=greenwhich_argument(astro_angles[m],wave);
      if(nodal==1) {
        Vu=nodal_phase(astro_angles[m],wave);
        f=nodal_factor(astro_angles[m],wave.formula);
        }
      else {
        Vu=0.;
        f=1.;
        }
      double V=V0+Vu;
       
      T cs=cos(V);
      T sn=sin(V);
      serie[m]+=f*(real(constants[k])*cs-imag(constants[k])*sn);
      }
    }
    
  delete[] astro_angles;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_predictionTS(float *serie, double *time, int n, spectrum_t s, fcomplex *constants, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=harmonic_predictionTS_template(serie, time, n,  s, constants, nodal);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_predictionTS(double *serie, double *time, int n, spectrum_t s, fcomplex *constants, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  
  status=harmonic_predictionTS_template(serie, time, n,  s, constants, nodal);
  
  return(status);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_predictionTS(double *serie, double *time, int n, mgr_data_t *data, int ndata, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  fcomplex *constants=new fcomplex[ndata];
  spectrum_t s(ndata);
  
  s.n=ndata;
  for(int k=0;k<s.n;k++) {
    constants[k]=polar(data[k].amp,-data[k].phi*M_PI/180.);
    s.waves[k]=data[k].constituent;
    }
  status=harmonic_predictionTS(serie, time, n, s, constants, nodal);
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_predictionTS(double *serie, double *time, int n, hconstant_t & data, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  spectrum_t s=*data.s;
  fcomplex *constants=new fcomplex[s.n];
  
  for(int k=0;k<s.n;k++) {
    constants[k]=polar((double) data.a[k],-data.G[k]*M_PI/180.);
    }
  status=harmonic_predictionTS(serie, time, n, s, constants, nodal);
  
  delete[] constants;
  
  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_predictionTS(tseries_t & serie, date_t start, date_t final, double dt, hconstant_t & mgr, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  int status;
  double *time;
  const spectrum_t s=*mgr.s;
  fcomplex *constants=new fcomplex[s.n];
  
  double startd=cnes_time(start,'d');
  double elapsed=time_elapsed(start,final);
  
  int nframes=1+(int) ceil(elapsed/dt);
  
  time=new double[nframes];
  const double incr=dt/24./3600.;
  
//   for(int k=0;k<nframes;k++) {
//     time[k]=startd+(double) k*incr;
//     }
  time[0]=startd;
  for(int k=1;k<nframes;k++) {
    time[k]=time[k-1]+incr;
    }
    
  for(int k=0;k<s.n;k++) {
    constants[k]=polar((double) mgr.a[k],-mgr.G[k]*M_PI/180.);
    }
    
  serie.n=nframes;
  serie.nparam=1;
  serie.x=new double*[1];
  serie.x[0]=new double[serie.n];
  serie.t=time;
    
  status=harmonic_predictionTS(serie.x[0], time, serie.n, s, constants, nodal);
  
//   delete[] time;
  delete[] constants;

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_prediction(double *buffer,double time, int n, spectrum_t s, fcomplex **constants, int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/* *----------------------------------------------------------------------

Purpose : Perform a tidal prediction at a given time for different positions

time is CNES time (origin 1950/01/01 00:00:00);
time units are decimal days;

----------------------------------------------------------------------*/
{
  int k,m;
  tidal_wave wave;
  double V,V0;
  double f,Vu;
  date_t start;
  extern date_t cnesdate(double t, char fmt);
  float cs,sn;
  astro_angles_t astro_angles;

  for(m=0; m<n; m++) {
    buffer[m]=0.0;
    }

  start=poctime_getdatecnes(time, 'd');
  init_argument(&astro_angles,start);

  for(k=0; k < s.n; k++) {
    wave=s.waves[k];
    V0=greenwhich_argument(astro_angles,wave);
    if(nodal==1) {
      Vu=nodal_phase(astro_angles,wave);
      f=nodal_factor(astro_angles,wave.formula);
      }
    else {
      Vu=0.;
      f=1.;
      }
/*     V=V0+Vu-G[k][m]; */
    V=V0+Vu;
    cs=cos(V);
    sn=sin(V);
    for(m=0; m<n; m++) {
/*       buffer[m]+=f*a[k][m]*cos(V); */
/*
      buffer[m]+=f*(constants[k][m].r*cs-constants[k][m].i*sn);
*/
      buffer[m]+=f*(real(constants[k][m])*cs-imag(constants[k][m])*sn);
      }
    }

  return(0);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  int harmonic_prediction(const astro_angles_t &astro_angles,double *buffer,double time,int n,const spectrum_t &s,hconstant_t *constants,int nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///Perform a tidal prediction at a given time for different points
/**
\param *buffer of n points
\param time in seconds, referenced to the last call of init_argument()
\param n number of points
\param s list of waves
\param[in,out] *constants array of n hconstants_t, initialised for s.n waves. Will be corrected if it mixes polar and complex elements.
\param nodal whether nodal corrections are carried out. Default:1
\returns 0 on success and EINVAL on argument mismatch
*/
/*----------------------------------------------------------------------------*/
{
  int k,m;//wave index,point index
  tidal_wave wave;
  double V,V0;
  double f,Vu;
  date_t start;
  float cs,sn;
  int ispolar,iscomplex,allPolar=1,allComplex=1;
  
/*------------------------------------------------------------------------------
  Ensure compatibility for mixed polar and complex constants */
  
  for(m=0; m<n; m++) {
    /* for all points */
    
    /* also reset summation of wave contribution */
    buffer[m]=0.0;
    
    if(constants[m].size!=s.n)
      TRAP_ERR_RETURN(EINVAL,1,"There are %d constants for point number %d but %d waves",constants[m].size,m,s.n);
    
    ispolar   = (constants[m].a!=NULL) && (constants[m].G!=NULL);
    iscomplex =  constants[m].z!=NULL ;
    
    if(!ispolar && !iscomplex){
      buffer[m]=NAN;
      TRAP_ERR_RETURN(EINVAL,1,"No constants for point number %d",m);
      }
    
    if(!ispolar)
      allPolar=0;
    
    if(!iscomplex)
      allComplex=0;
    
    }
  
  if(allPolar==0 && allComplex==0){
    /* if constants mixes polar and complex elements */
    
#if 1
    
    /* set all complex to polar */
    for(m=0; m<n; m++) {
      if(constants[m].a==NULL || constants[m].G==NULL)
        constants[m].set_polar(1);
      }
    
    allPolar=1;
    
#else
    
    /* set all polar to complex (for testing) */
    for(m=0; m<n; m++) {
      if(constants[m].z==NULL)
        constants[m].set_complex();
      }
    
    allComplex=1;
    
#endif
    
    }
  
/*------------------------------------------------------------------------------
  PREDICTION */
  
  /* compute astronomic angles */
  astro_angles_t shifted_astro_angles;
  astronomic_angle(&shifted_astro_angles,astro_angles.t_julian+time/(d2s*jc2d),0);
  
  FILE *F=0;
//   F=fopen((string(__func__)+".dat").c_str(),"a");
  if(F!=0) fprintf(F,"%g",shifted_astro_angles.t_julian);
  
  for(k=0; k < s.n; k++) {
    /* For all waves */
    
    /* compute greenwhich argument, eventually with nodal correction */
    wave=s.waves[k];
    bool waveIsM0=(strcmp(wave.name,"M0")==0);
    
    V0=greenwhich_argument(shifted_astro_angles,wave);
    
    if(nodal==1) {
      Vu=nodal_phase(shifted_astro_angles,wave);
      f=nodal_factor(shifted_astro_angles,wave.formula);
      }
    else {
      Vu=0.;
      f=1.;
      }
    
    V=V0+Vu;
    if(F!=0) fprintf(F," %g %g",f,V);
    
    /* add wave contribution to all the points */
    
    if(allPolar){
      
      #pragma omp parallel for if(n>10000)
      for(m=0; m<n; m++) { 
        buffer[m]+=f*constants[m].a[k]*cos(V-constants[m].G[k]*d2r);
        if(waveIsM0)
          buffer[m]-=constants[m].a[k]*cos(-constants[m].G[k]*d2r);
        }
      
      }
    else{
      cs=cos(V);
      sn=sin(V);
      
      #pragma omp parallel if(n>10000)
      {
      complex<float> z;
      
      #pragma omp for
      for(m=0; m<n; m++) {
        z=constants[m].z[k];
        buffer[m]+=f*(real(z)*cs-imag(z)*sn);
        if(waveIsM0)
          buffer[m]-=real(z);
        }
      
      }//EO omp parallel
      
      }
    
    }
  
  if(F!=0){
    fprintf(F,"\n");
    fclose(F);
    }
  
  return 0;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

template<typename T> void harmonic_prediction_fortran_template(const astro_angles_t *astro_angles,double *buffer,const double *time,int *n,const spectrum_t *s,const T *amplitudes,const T *degrees,int *nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*----------------------------------------------------------------------------*/
///template Fortran wrapper for harmonic_prediction(const astro_angles_t &astro_angles, double *buffer,double time,int n,const spectrum_t &s, hconstant_t *constants, int nodal)
/**
\param *buffer of n points
\param *time in seconds, referenced to the last call of init_argument()
\param *n number of points
\param *s list of waves
\param *amplitudes array of (s.n,n) amplitudes
\param *degrees array of (s.n,n) phase delays
\param *nodal whether nodal corrections are carried out.
*/
/*----------------------------------------------------------------------------*/
{
  int i,j;//point and wave indexes
  hconstant_t *constants;
  
  /** \bug : this uses hconstant_t, which in turn refers to a float */
  constants=new hconstant_t[*n];
  
  for(i=0;i<*n;i++){//for all points
    constants[i].init_polar(s->n);
    for(j=0;j<s->n;j++){//for all waves
      constants[i].a[j]=amplitudes[j* *n+i];
      constants[i].G[j]=degrees[j* *n+i];
      }
    }
  
  harmonic_prediction(*astro_angles,buffer,*time,*n,*s,constants,*nodal);
  
  for(i=0;i<*n;i++)
    constants[i].destroy();
  delete[]constants;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void harmonic_prediction4_(const astro_angles_t *astro_angles,double *buffer,const double *time,int *n,const spectrum_t *s,int *nwaves,const float *amplitudes,const float *degrees,int *nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  harmonic_prediction_fortran_template(astro_angles,buffer,time,n,s,amplitudes,degrees,nodal);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

extern "C" void harmonic_prediction8_(const astro_angles_t *astro_angles,double *buffer,const double *time,int *n,const spectrum_t *s,int *nwaves,const double *amplitudes,const double *degrees,int *nodal)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  harmonic_prediction_fortran_template(astro_angles,buffer,time,n,s,amplitudes,degrees,nodal);
}
