
/**************************************************************************

  T-UGOm hydrodynamic ocean model, 2006-2012

  Unstructured Ocean Grid initiative

Contributors:

  Florent Lyard      LEGOS/CNRS, Toulouse, France
  Cyril Nguyen       LA/CNRS,    Toulouse, France
  Laurent Roblou     LEGOS/CNRS, Toulouse, France
  Damien Allain      LEGOS/CNRS, Toulouse, France
  
  David Greenberg    Bedford Institute of Oceanography, Halifax, Canada
  Frédéric Dupont    Université de Laval à Québec, Canada
  
  Yoann Le Bars      PhD, LEGOS, Toulouse, France
  Yves Soufflet      Post-doctorant, LEGOS, Toulouse, France
  Clement Mayet      PhD, LEGOS, Toulouse, France

***************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <cmath>

#define MAX(a,b)  (((a) > (b)) ? (a) : (b))

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double svan(double s4, double t4, double p04, double *sigma)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  /* Initialized data */

  static double r3500 = 1028.1063;
  static double r4    = 4.8314e-4;
  static double dr350 = 28.106331;

  /* System generated locals */
  double ret_val;

  double r1, r2, r3;
  double e4, b4, bw;
  double d4, bb1, a4, aw, aa1;
  double k0, kk, kw;

  /* Local variables */
  double dvan, dr35p;
  double c4;
  double p4, dk, k35;
  double pk;
  double sr;
  double gam, sig, sva, v350p;

/**----------------------------------------------------------------------

 Specific volume anomaly (steric anomaly) based on 1980 equation
 of state for seawater and 1978 practerical salinity scale.
 References:
 Millero, et al (1980) Deep-Sea Res.,27A,255-264.
 Millero and Poisson 1981,Deep-Sea Res.,28A PP 625-629.
 Both above references are also found in UNESCO REPORT 38 (1981)

 UNITS:
     Pressure	     P04    Decibars
     Temperature     T4     Deg Celsius (IPTS-68)
     Salinity	     S4     (IPSS-78)
     Spec. Vol. Ana. SVAN   m**3/kg *1.0E-8
     Density Ana.    SIGMA  kg/m**3

 Check value: SVAN=981.3021 E-8 m**3/kg. For S = 40 (IPSS-78),
 T = 40 Deg C, P0= 10000 Decibars.
 Check value: SIGMA = 59.82037  kg/m**3. For S = 40 (IPSS-78),
 T = 40 Deg C, P0= 10000 Decibars.

 ********************
 DATA

   RR4 is refered to as C in notation of Millero and Poisson 1981
   Convert pressure to bars and take square root salinity.

----------------------------------------------------------------------*/

  p4 = p04 / 10.f;
  sr = sqrt(fabs(s4));

/*----------------------------------------------------------------------
   Pure water density at atmospheric pressure
   BIGG P.H.,(1967) BR. J. APPLIED PHYSICS 8 PP 521-537. */

  r1 = ((((t4 * 6.536332e-9f - 1.120083e-6f) * t4 + 1.001685e-4f) * t4 -
         .00909529f) * t4 + .06793952f) * t4 - 28.263737f;

/*----------------------------------------------------------------------
  seawater density atm press.
  coefficients involving salinity
  RR2 = A   in notation of Millero and Poisson 1981
-----------------------------------------------------------------------*|
  r2 = (((t4 * 5.3875e-9f - 8.2467e-7f) * t4 + 7.6438e-5f) * t4 -
      .0040899f) * t4 + .824493f;

/*----------------------------------------------------------------------
  RR3 = B4  in notation of Millero and Poisson 1981 */
  r3 = (t4 * -1.6546e-6f + 1.0227e-4f) * t4 - .00572466f;

/*----------------------------------------------------------------------
  International one-atmosphere equation of state of seawater */
  sig = (r4 * s4 + r3 * sr + r2) * s4 + r1;

/*----------------------------------------------------------------------
 Specific volume at atmospheric pressure */
  v350p = 1.f / r3500;
  sva = -sig * v350p / (r3500 + sig);

  *sigma = sig + dr350;

/*----------------------------------------------------------------------
  Scale specific vol. anamoly to normally reported units */
  ret_val = sva * 1e8f;

  if(p4 == 0.f) {
    return (ret_val);
  }

/*----------------------------------------------------------------------
  New High Pressure Equation of Sate for Seawater
                              
  Millero, el al., 1980 DSR 27A, pp 255-264
  Constant notation follows article
-----------------------------------------------------------------------*|

/*----------------------------------------------------------------------
 Compute compression terms */
  e4 = (t4 * 9.1697e-10f + 2.0816e-8f) * t4 - 9.9348e-7f;
  bw = (t4 * 5.2787e-8f - 6.12293e-6f) * t4 + 3.47718e-5f;
  b4 = bw + e4 * s4;

  d4 = 1.91075e-4f;
  c4 = (t4 * -1.6078e-6f - 1.0981e-5f) * t4 + .0022838f;
  aw = ((t4 * -5.77905e-7f + 1.16092e-4f) * t4 + .00143713f) * t4 - .1194975f;
  a4 = (d4 * sr + c4) * s4 + aw;

  bb1 = (t4 * -5.3009e-4f + .016483f) * t4 + .07944f;
  aa1 = ((t4 * -6.167e-5f + .0109987f) * t4 - .603459f) * t4 + 54.6746f;
  kw = (((t4 * -5.155288e-5f + .01360477f) * t4 - 2.327105f) * t4 +
        148.4206f) * t4 - 1930.06f;
  k0 = (bb1 * sr + aa1) * s4 + kw;

/*----------------------------------------------------------------------
  Evaluate pressure polynomial
  -----------------------------------------------------
  K Equals the secant bulk modulus of seawater
  DK=K(S,T,P)-K(35,0,P)
  K35=K(35,0,P)
----------------------------------------------------------------------*/
  dk = (b4 * p4 + a4) * p4 + k0;
  k35 = (p4 * 5.03217e-5f + 3.359406f) * p4 + 21582.27f;
  gam = p4 / k35;
  pk = 1.f - gam;
  sva = sva * pk + (v350p + sva) * p4 * dk / (k35 * (k35 + dk));

/*----------------------------------------------------------------------
  SCALE SPECIFIC VOL. ANAMOLY TO NORMALLY REPORTED UNITS */
  ret_val = sva * 1e8f;
  v350p *= pk;

/*----------------------------------------------------------------------
 Compute density anamoly with respect to 1000.0 kg/m**3
  1) DR350: Density anamoly at 35 (IPSS-78),
                 0 Deg. C and 0 decibars
  2) DR35P: Density anamoly at 35 (IPSS-78),
                 0 Deg. C, Pres. variation
  3) DVAN : Density anamoly variations involving specific
      volume anamoly		
                              
 Check values: SIGMA = 59.82037 kg/m**3 	
 for S = 40 (IPSS-78), T = 40 Deg C, P0= 10000 Decibars.
----------------------------------------------------------------------*/
  dr35p = gam / v350p;
  dvan = sva / (v350p * (v350p + sva));
  *sigma = dr350 + dr35p - dvan;
  return ret_val;
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double physic_theta(double s4, double t04, double p04, double pr)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  /* System generated locals */
  double ret_val;

  /* Local variables */
  double h4, p4, q4, t4, xk;
  extern double atg(double, double, double);

/**----------------------------------------------------------------------

 To compute local potential temperature at PR using	
 Bryden 1973 Polynomial for adiabatic lapse rate and	
 Runge-Kutta 4th order integration algorithm.
 Ref: Bryden,H.,1973,Deep-Sea Res.,20,401-408;
 Fofonoff,N.,1977,Deep-Sea Res.,24,489-491

 Units:
     Pressure         P04     Decibars
     Temperature      T04     Deg Celsius (IPTS-68)
     Salinity         S4      (IPSS-78)
     Reference PRS    PR      Decibars
     Potential TMP.   THETA   Deg Celsius

 Checkvalue:
       THETA= 36.89073 C,S=40 (IPSS-78),
       T0=40 DEG C,P0=10000 Decibars,PR=0 decibars

 Set up intermediate temperature and pressure variables.

----------------------------------------------------------------------*/

  p4 = p04;
  t4 = t04;
  h4 = pr - p4;
  xk = h4 * atg(s4, t4, p4);

  t4 += xk * .5f;
  q4 = xk;
  p4 += h4 * .5f;
  xk = h4 * atg(s4, t4, p4);

  t4 += (xk - q4) * .29289322f;
  q4 = xk * .58578644f + q4 * .121320344f;
  xk = h4 * atg(s4, t4, p4);

  t4 += (xk - q4) * 1.707106781f;
  q4 = xk * 3.414213562f - q4 * 4.121320344f;
  p4 += h4 * .5f;
  xk = h4 * atg(s4, t4, p4);

  ret_val = t4 + (xk - q4 * 2.f) / 6.f;

  return (ret_val);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double atg(double s4, double t4, double p4)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double ret_val;
  double ds;

/*----------------------------------------------------------------------------*/
/**
 Adiabatic temperature gradient Deg c per decibar
 Ref: Bryden, H., 1973,Deep-Sea Res.,20,401-408,
      http://dx.doi.org/10.1016/0011-7471(73)90063-6

 Units:
     Pressure      P4     Decibars
     Temperature   T4     Deg Celsius(IPTS-68)
     Salinity      S4     (IPSS-78)
     Adiabatic     ATG    Deg. C/Decibar

 Checkvalue:

 ATG=3.255976E-4 C/DBAR FOR S=40 (IPSS-78),T=40 Deg C,P0=10000 Decibars

*/
/*----------------------------------------------------------------------------*/

  ds = s4 - 35.f;
  ret_val =
        (((t4 * -2.1687e-16f + 1.8676e-14f) * t4 - 4.6206e-13f) * p4 +
         ((t4 * 2.7759e-12f  - 1.1351e-10f) * ds +
         ((t4 * -5.4481e-14f + 8.733e-12f)  * t4 - 6.7795e-10f) * t4 +
          1.8741e-8f)) * p4 + (t4 * -4.2393e-8f + 1.8932e-6f) * ds +
        ((t4 * 6.6228e-10f - 6.836e-8f) * t4 + 8.5258e-6f) * t4 + 3.5803e-5f;

  return (ret_val);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double water_density_00(double tf, double sf)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  double ret_val, d__1;
  double rhof;

/*----------------------------------------------------------------------

  version(03/02/90), origin ???
  this function computes density

----------------------------------------------------------------------*/

/*----------------------------------------------------------------------
  salinity dependance*/
  rhof = sf * sf * sf * 6.76786136e-6 - sf * sf * 4.8249614e-4 +
         sf * 0.814876577 - 0.22584586;

/*----------------------------------------------------------------------
  temperature rectification*/
  rhof *= tf * tf * tf * 1.667e-8 - tf * tf * 8.164e-7 + tf * 1.803e-5;

/*----------------------------------------------------------------------
  temperature dependance*/
  rhof = rhof + 1. - tf * tf * tf * 1.0843e-6 + tf * tf * 9.8185e-5 -
         tf * .004786;
/*----------------------------------------------------------------------
  salinity rectification*/
  rhof *= sf * sf * sf * 6.76786136e-6 - sf * sf * 4.8249614e-4 +
          sf * .814876577 + .03895414;

/*----------------------------------------------------------------------
  Computing 2nd power */
  d__1 = tf - 3.98;
  rhof -= d__1 * d__1 * (tf + 283.) / ((tf + 67.26) * 503.57);

/*----------------------------------------------------------------------
  modified to return full density (not anomaly)*/
  ret_val = rhof +1000.;

  return (ret_val);
}



/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double water_density_01(double t, double s, double p)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
/*------------------------------------------------------------------------------

  Algorithms for Density, Potential Temperature, Conservative Temperature, and the
  Freezing Temperature of Seawater
  DAVID R. JACKETT et al., JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY, 2006
  http://journals.ametsoc.org/doi/abs/10.1175/JTECH1946.1

  An Assessment of Orthobaric Density in the Global Ocean. 
  Trevor J. McDougall and David R. JackettJ. Phys. Oceanogr., 35, 2054–2075.
  doi: http://dx.doi.org/10.1175/JPO2796.1

  A Neutral Density Variable for the World’s Oceans
  David R. Jackett and Trevor J. McDougall, JOURNAL OF PHYSICAL OCEANOGRAPHY, 1997
  http://journals.ametsoc.org/doi/abs/10.1175/1520-0485%281997%29027%3C0237:ANDVFT%3E2.0.CO;2
  
------------------------------------------------------------------------------*/

{
double    c1_jmfwg                           \
         ,c2_jmfwg                           \
         ,c3_jmfwg                           \
         ,c4_jmfwg                           \
         ,c5_jmfwg                           \
         ,c6_jmfwg                           \
         ,c7_jmfwg                           \
         ,c8_jmfwg                           \
         ,c9_jmfwg                           \
         ,c10_jmfwg                          \
         ,c11_jmfwg                          \
         ,c12_jmfwg                          \
         ,c14_jmfwg                          \
         ,c15_jmfwg                          \
         ,c16_jmfwg                          \
         ,c17_jmfwg                          \
         ,c18_jmfwg                          \
         ,c19_jmfwg                          \
         ,c20_jmfwg                          \
         ,c21_jmfwg                          \
         ,c22_jmfwg                          \
         ,c23_jmfwg                          \
         ,c24_jmfwg                          \
         ,c25_jmfwg;
          
  double rho, Pn,Pd,s32,t3;
  
  s=MAX(s,0.0);
  s32=exp(3./2.0*log(s));
  t3=t*t*t;
  
/**--------------------------------------------------------------------------------
  Coefficients of the Jackett et al. (2006) equation of state */
  c1_jmfwg=  9.9984085444849347e+2;
  c2_jmfwg=  7.3471625860981584e+0;
  c3_jmfwg= -5.3211231792841769e-2;
  c4_jmfwg=  3.6492439109814549e-4;
  c5_jmfwg=  2.5880571023991390e+0 ;
  c6_jmfwg= -6.7168282786692355e-3;
  c7_jmfwg=  1.9203202055760151e-3;
  c8_jmfwg=  1.1798263740430364e-2;
  c9_jmfwg=  9.8920219266399117e-8;
  c10_jmfwg= 4.6996642771754730e-6;
  c11_jmfwg=-2.5862187075154352e-8;
  c12_jmfwg=-3.2921414007960662e-12;
  c14_jmfwg= 7.2815210113327091e-3;
  c15_jmfwg=-4.4787265461983921e-5;
  c16_jmfwg= 3.3851002965802430e-7;
  c17_jmfwg= 1.3651202389758572e-10;
  c18_jmfwg= 1.7632126669040377e-3;
  c19_jmfwg=-8.8066583251206474e-6;
  c20_jmfwg=-1.8832689434804897e-10;
  c21_jmfwg= 5.7463776745432097e-6;
  c22_jmfwg= 1.4716275472242334e-9;
  c23_jmfwg= 6.7103246285651894e-6;
  c24_jmfwg=-2.4461698007024582e-17;
  c25_jmfwg=-9.1534417604289062e-18;


  
  Pn= c1_jmfwg + t*(c2_jmfwg+t*(c3_jmfwg+t*c4_jmfwg)) \
               + s*(c5_jmfwg+t*c6_jmfwg+s*c7_jmfwg)   \
               + p*(c8_jmfwg+t*t*c9_jmfwg+s*c10_jmfwg+p*(c11_jmfwg+t*t*c12_jmfwg));

  Pd= 1.0 + t*(c14_jmfwg+t*(c15_jmfwg+t*(c16_jmfwg+t*c17_jmfwg)))  \
          + s*(c18_jmfwg+t*c19_jmfwg+t*t*t*c20_jmfwg)                 \
          + sqrt(s*s*s)*(c21_jmfwg+t*t*c22_jmfwg)                           \
          + p*(c23_jmfwg+p*(t*t*t*c24_jmfwg+p*t*c25_jmfwg));

  rho=Pn/Pd;
  
  return(rho);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double water_density_02(double t, double s)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
double    c1_jmfwg                           \
         ,c2_jmfwg                           \
         ,c3_jmfwg                           \
         ,c4_jmfwg                           \
         ,c5_jmfwg                           \
         ,c6_jmfwg                           \
         ,c7_jmfwg                           \
         ,c8_jmfwg                           \
         ,c9_jmfwg                           \
         ,c10_jmfwg                          \
         ,c11_jmfwg                          \
         ,c12_jmfwg                          \
         ,c13_jmfwg                          \
         ,c14_jmfwg                          \
         ,c15_jmfwg                          \
         ,c16_jmfwg;
          
  double rho, Pn, Pd;
  
/**-----------------------------------------------------------------------------
  Coefficients of the Jackett et al. (2005) equation of state */
  c1_jmfwg=  1.0023063688892480e+3;
  c2_jmfwg=  2.2280832068441331e-1;
  c3_jmfwg=  8.1157118782170051e-2;
  c4_jmfwg= -4.3159255086706703e-4;
  c5_jmfwg= -1.0304537539692924e-4;
  c6_jmfwg= -3.1710675488863952e-3;
  c7_jmfwg= -1.7052298331414675e-7;
  
  c8_jmfwg=  4.3907692647825900e-5;
  c9_jmfwg=  7.8717799560577725e-5;
  c10_jmfwg=-1.6212552470310961e-7;
  c11_jmfwg=-2.3850178558212048e-9;
  c12_jmfwg=-5.1268124398160734e-4;
  c13_jmfwg= 6.0399864718597388e-6;
  c14_jmfwg=-2.2744455733317707e-9;
  c15_jmfwg=-3.6138532339703262e-5;
  c16_jmfwg=-1.3409379420216683e-9;

  s=MAX(s,0.0);
  
  Pn= c1_jmfwg + t*(c2_jmfwg+t*(c3_jmfwg+t*c4_jmfwg)) \
               + s*(c5_jmfwg+t*c6_jmfwg+s*c7_jmfwg);

  Pd= 1.0 + t*(c8_jmfwg+t*(c9_jmfwg+t*(c10_jmfwg+t*c11_jmfwg)))  \
          + s*(c12_jmfwg+t*c13_jmfwg+t*t*t*c14_jmfwg)                 \
          + sqrt(s*s*s)*(c15_jmfwg+t*t*c16_jmfwg);

  rho=Pn/Pd;
  
  return(rho);
}

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double water_density(double t, double s, double z, int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
  extern double sigma0 ( float **ptem, float **psal, int kpi, int  kpj)  ;
  double p,rho;
  float **ptem, **psal;
  int kpi, kpj;
 
  switch(option) {
    case 0:
      rho=water_density_00(t, s);
      break;
      
    case 1:
/*------------------------------------------------------------------------------
      t potential temperature (celsius), s salinity (psu), p pressure (dbar) */
      p=(101300+z*9.81*1025)*0.0001;
// i.e. p=(101300+z*9.81*1025)*0.00001*10;
      rho=water_density_01(t, s, p);
      break;
      
     case 2:
      rho=water_density_02(t, s);
      break;
    }
    
//  rho=sigma0 ( ptem, psal, kpi,  kpj)  ;
  return(rho);
}


/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/

  double water_density_check(int option)

/*XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX*/
{
/*------------------------------------------------------------------------------

  Algorithms for Density, Potential Temperature, Conservative Temperature, and the
  Freezing Temperature of Seawater
  DAVID R. JACKETT et al., JOURNAL OF ATMOSPHERIC AND OCEANIC TECHNOLOGY, 2006
  
------------------------------------------------------------------------------*/

  double p,T,S,rho,r;
  
  S=35; T=25; p=2000; rho=1031.65056056576;
  r=water_density_01(T, S, p);
  printf("S=%lf T=%lf p=%lf\n",S,T,p);
  printf("r=%15lf ref=%15lf chk=%15lf\n",r,rho,r-rho);
  
  S=20; T=20; p=1000; rho=1017.72886801964;
  r=water_density_01(T, S, p);
  printf("S=%lf T=%lf p=%lf\n",S,T,p);
  printf("r=%15lf ref=%15lf chk=%15lf\n",r,rho,r-rho);
  
  S=40; T=12; p=8000; rho=1062.95279820631;
  r=water_density_01(T, S, p);
  printf("S=%lf T=%lf p=%lf\n",S,T,p);
  printf("r=%15lf ref=%15lf chk=%15lf\n",r,rho,r-rho);
  
  return(rho);
  
}
